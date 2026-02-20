from __future__ import annotations

from typing import Dict, List, Tuple, Optional

import polars as pl
import numpy as np

from ..utils.logging import log_info
from .frame_crosstalk import estimate_crosstalk_matrix, invert_and_correct
from .junction_model import estimate_junction_expectation, junction_llr


def _profiles_to_tran_df(prof: pl.DataFrame, tran_id: Optional[str] = None) -> pl.DataFrame:
    df = prof
    if tran_id is not None:
        df = df.filter(pl.col('tran_id') == tran_id)
    return df.rename({'pos': 'tran_start', 'count': 'counts'})


def _frame_counts(df: pl.DataFrame) -> np.ndarray:
    if df.is_empty():
        return np.zeros((0, 3), dtype=float)
    # aggregate per-codon position: use pos modulo 3
    f = (
        df.with_columns((pl.col('tran_start') % 3).alias('frame'))
          .group_by('frame')
          .agg(pl.col('counts').sum().alias('sum'))
          .sort('frame')
    )
    out = np.zeros((1, 3), dtype=float)
    for i in range(min(3, f.height)):
        out[0, int(f['frame'][i])] = float(f['sum'][i])
    return out


def compute_exon_chunk_metrics(
    profiles: pl.DataFrame,
    feature_row: dict,
    fmap: pl.DataFrame,
) -> dict:
    """Compute metrics for an exon_chunk by aggregating over transcripts that use it.

    We use transcript ranges from feature_map to translate genomic chunk to transcript coordinates per transcript.
    """
    fid = feature_row['feature_id']
    # Pull all transcript mappings that include this feature
    rows = fmap.filter(pl.col('feature_chain').list.contains(fid))
    cov_total = 0.0
    cov_f0 = 0.0
    nzc = 0.0
    n_codons = 0
    for r in rows.iter_rows(named=True):
        tx = r['transcript_id']
        chain: List[str] = r['feature_chain']
        i = chain.index(fid)
        tr_starts: List[int | None] = r['tran_ranges_start']
        tr_ends: List[int | None] = r['tran_ranges_end']
        ts, te = tr_starts[i], tr_ends[i]
        if ts is None or te is None:
            continue
        tx_df = _profiles_to_tran_df(profiles, tran_id=tx)
        seg = tx_df.filter((pl.col('tran_start') >= ts) & (pl.col('tran_start') < te))
        if seg.is_empty():
            continue
        fcounts = _frame_counts(seg)
        cov_total += float(seg['counts'].sum())
        cov_f0 += float(fcounts[0, 0])
        # NZC: fraction codons with nonzero in any frame within this segment (approx by positions/3)
        npos = te - ts
        n_cod = max(1, npos // 3)
        nzc += min(1.0, float(seg.filter(pl.col('counts') > 0).height) / float(npos or 1))
        n_codons += 1
    hrf = (cov_f0 / max(1.0, cov_total - cov_f0)) if cov_total > 0 else 0.0
    return {
        'cov_total_mean': cov_total,
        'cov_f0_mean': cov_f0,
        'nzc': (nzc / n_codons) if n_codons > 0 else 0.0,
        'hrf': hrf,
        # periodicity_p placeholder: compute elsewhere if needed
        'periodicity_p': 1.0,
    }


def compute_tis_tts_metrics(
    profiles: pl.DataFrame,
    feature_row: dict,
    fmap: pl.DataFrame,
    sru_range: int = 15,
) -> dict:
    from ..core.scoring import sru_score
    fid = feature_row['feature_id']
    rows = fmap.filter(pl.col('feature_chain').list.contains(fid))
    up_scores = []
    down_scores = []
    for r in rows.iter_rows(named=True):
        tx = r['transcript_id']
        chain: List[str] = r['feature_chain']
        i = chain.index(fid)
        pos_list: List[int | None] = r.get('tran_pos') or []
        if i >= len(pos_list):
            continue
        tpos = pos_list[i]
        if tpos is None:
            continue
        tx_df = _profiles_to_tran_df(profiles, tran_id=tx)
        # Convert to expected columns for sru_score
        # Already has tran_start, counts
        up = sru_score(int(tpos), tx_df, sru_range, 0)
        down = sru_score(int(tpos), tx_df, sru_range, 1)
        up_scores.append(up)
        down_scores.append(down)
    return {
        'sru_up': float(np.mean(up_scores)) if up_scores else 0.0,
        'ramp_symmetry': float(np.mean(up_scores) - np.mean(down_scores)) if up_scores and down_scores else 0.0,
        'sru_down': float(np.mean(down_scores)) if down_scores else 0.0,
        'step_symmetry': float(np.mean(down_scores) - np.mean(up_scores)) if up_scores and down_scores else 0.0,
        'local_periodicity': 0.0,
    }


def compute_junction_metrics(
    junction_features: pl.DataFrame,
    genome_bam_splits: pl.DataFrame,
    expectation_model: pl.DataFrame,
) -> pl.DataFrame:
    # For now, single global bin
    if expectation_model.is_empty():
        exp_frac = 0.0
        var = 1e-6
    else:
        exp_frac = float(expectation_model['E_split_fraction'][0])
        var = float(expectation_model['var'][0])
    # Aggregate observed splits per junction key
    g = genome_bam_splits.group_by(['chr', 'donor_pos', 'acceptor_pos', 'strand']).agg(pl.col('count').sum().alias('obs_split'))
    # Join to features
    jf = junction_features.join(g, on=['chr', 'donor_pos', 'acceptor_pos', 'strand'], how='left').fill_null(0)
    # Estimate total reads near junction (placeholder): use obs_split as proxy total
    jf = jf.with_columns(pl.col('obs_split').alias('total_reads'))
    # LLR
    jf = jf.with_columns(
        pl.struct(['obs_split', 'total_reads']).map_elements(lambda s: junction_llr(int(s['obs_split']), exp_frac, int(s['total_reads']), var)).alias('split_llr')
    )
    return jf


def feature_metrics(
    profiles_parquet: str,
    feature_parquet: str,
    feature_map_parquet: str,
    genome_bam_splits: Optional[pl.DataFrame] = None,
    out_parquet: str = 'feature_metrics.parquet',
    sru_range: int = 15,
) -> str:
    """Compute per-feature metrics and write to Parquet.

    genome_bam_splits: optional DataFrame containing split junction counts with columns
      chr, donor_pos, acceptor_pos, strand, count
    """
    profiles = pl.read_parquet(profiles_parquet)
    feats = pl.read_parquet(feature_parquet)
    fmap = pl.read_parquet(feature_map_parquet)

    # Crosstalk correction stub: left as a no-op here; can be extended to per-length
    # Expectation model stub
    exp_model = pl.DataFrame({'bin_id': [0], 'E_split_fraction': [0.0], 'var': [1e-6]})

    rows: List[dict] = []
    # Split features by type
    chunks = feats.filter(pl.col('feature_type') == 'exon_chunk')
    juncs = feats.filter(pl.col('feature_type') == 'junction')
    tis = feats.filter(pl.col('feature_type') == 'TIS')
    tts = feats.filter(pl.col('feature_type') == 'TTS')

    # Exon chunks
    for r in chunks.iter_rows(named=True):
        m = compute_exon_chunk_metrics(profiles, r, fmap)
        rows.append({**r, **m})

    # TIS/TTS
    for r in tis.iter_rows(named=True):
        m = compute_tis_tts_metrics(profiles, r, fmap, sru_range=sru_range)
        rows.append({**r, **m})
    for r in tts.iter_rows(named=True):
        m = compute_tis_tts_metrics(profiles, r, fmap, sru_range=sru_range)
        rows.append({**r, **m})

    # Junctions
    if genome_bam_splits is not None and not genome_bam_splits.is_empty():
        jm = compute_junction_metrics(juncs, genome_bam_splits, exp_model)
        rows.extend([row for row in jm.iter_rows(named=True)])
    else:
        # No split info; add with zeros
        for r in juncs.iter_rows(named=True):
            rows.append({**r, 'split_llr': 0.0})

    out = pl.from_dicts(rows)
    out.write_parquet(out_parquet)
    log_info(f"Feature metrics written: {out_parquet}")
    return out_parquet
