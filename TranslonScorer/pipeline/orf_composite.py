from __future__ import annotations

from typing import Dict, List, Tuple

import polars as pl

DEFAULT_WEIGHTS = {
    'w_TIS': 1.0,
    'w_TTS': 1.0,
    'w_chunk': 0.5,
    'w_junc': 1.0,
    'w_per': 0.5,
    'w_cov': 0.25,
}


def _load_df(path: str) -> pl.DataFrame:
    if path.endswith('.parquet'):
        return pl.read_parquet(path)
    elif path.endswith('.csv') or path.endswith('.tsv'):
        sep = '\t' if path.endswith('.tsv') else ','
        return pl.read_csv(path, separator=sep)
    else:
        return pl.read_parquet(path)


def _features_for_orf(
    tran_id: str,
    orf_start: int,
    orf_stop: int,
    fmap_row: dict,
    feature_types: Dict[str, str],
) -> List[str]:
    """Return ordered list of feature_ids used by ORF along the transcript path.

    We include exon_chunk features whose transcript ranges overlap [orf_start, orf_stop),
    junction features that sit between included chunks, and TIS/TTS if their transcript
    positions match the ORF start/stop respectively.
    """
    chain: List[str] = fmap_row['feature_chain']
    tr_starts: List[int | None] = fmap_row['tran_ranges_start']
    tr_ends: List[int | None] = fmap_row['tran_ranges_end']
    tr_pos: List[int | None] = fmap_row.get('tran_pos') or [None] * len(chain)

    used: List[str] = []
    prev_in = False
    for i, fid in enumerate(chain):
        ftype = feature_types.get(fid, '')
        if ftype == 'exon_chunk':
            ts, te = tr_starts[i], tr_ends[i]
            if ts is None or te is None:
                prev_in = False
                continue
            if not (te <= orf_start or ts >= orf_stop):
                used.append(fid)
                prev_in = True
            else:
                prev_in = False
        elif ftype == 'junction':
            # Include junction if ORF spans across adjacent included chunks
            if prev_in:
                # if next chunk is also within ORF, include this junction
                # we approximate by checking neighbors' overlap
                nxt_in = False
                if i + 1 < len(chain):
                    nts, nte = (tr_starts[i + 1], tr_ends[i + 1])
                    if nts is not None and nte is not None:
                        nxt_in = not (nte <= orf_start or nts >= orf_stop)
                if nxt_in:
                    used.append(fid)
            prev_in = False
        elif ftype in ('TIS', 'TTS'):
            p = tr_pos[i]
            if p is not None:
                if ftype == 'TIS' and int(p) == int(orf_start):
                    used.append(fid)
                if ftype == 'TTS' and int(p) == int(orf_stop):
                    used.append(fid)
            prev_in = False
        else:
            prev_in = False
    return used


def orf_composite(
    orfs_path: str,
    feature_metrics_path: str,
    feature_map_path: str,
    out_parquet: str,
    weights: Dict[str, float] | None = None,
) -> str:
    """Aggregate per-feature metrics into composite ORF scores and write Parquet.

    orfs_path: Parquet/CSV with at least columns: orf_id, tran_id, start, stop, frame, type, length, (locus_id optional)
    feature_metrics_path: Parquet from feature-metrics
    feature_map_path: Parquet from features command
    """
    W = dict(DEFAULT_WEIGHTS)
    if weights:
        W.update(weights)

    orfs = _load_df(orfs_path)
    fmet = _load_df(feature_metrics_path)
    fmap = _load_df(feature_map_path)

    # Map feature_id -> metrics row, and feature_id -> type
    feat_type = {r[0]: r[1] for r in fmet.select(['feature_id', 'feature_type']).iter_rows()}
    frows = {r[0]: r for r in fmet.iter_rows(named=True)}

    # Build transcript map dict
    tmap = {}
    for r in fmap.iter_rows(named=True):
        tmap.setdefault(r['transcript_id'], []).append(r)

    out_rows: List[dict] = []
    for r in orfs.iter_rows(named=True):
        tx = r['tran_id']
        start = int(r['start'])
        stop = int(r['stop'])
        # find mapping row for this transcript (feature chain)
        maps = tmap.get(tx, [])
        if not maps:
            continue
        fmap_row = maps[0]
        fids = _features_for_orf(tx, start, stop, fmap_row, feat_type)
        # Aggregate
        s_TIS = s_TTS = per_term = cov_term = junc_term = chunk_term = 0.0
        for fid in fids:
            row = frows.get(fid)
            if not row:
                continue
            ftype = feat_type.get(fid, '')
            if ftype == 'TIS':
                s_TIS += float(row.get('sru_up', 0.0))
            elif ftype == 'TTS':
                s_TTS += float(row.get('sru_down', 0.0))
            elif ftype == 'junction':
                junc_term += float(row.get('split_llr', 0.0))
            elif ftype == 'exon_chunk':
                # coverage and periodicity proxies
                cov = float(row.get('cov_total_mean', 0.0))
                hrf = float(row.get('hrf', 0.0))
                cov_term += cov
                # per-term: transform HRF into bounded score
                per_term += (hrf / (1.0 + hrf)) if hrf >= 0 else 0.0
                chunk_term += 1.0
        score = (
            W['w_TIS'] * s_TIS +
            W['w_TTS'] * s_TTS +
            W['w_junc'] * junc_term +
            W['w_chunk'] * chunk_term +
            W['w_cov'] * cov_term +
            W['w_per'] * per_term
        )
        out_rows.append({**r, 'feature_ids': fids, 'score': score})

    out = pl.from_dicts(out_rows)
    out.write_parquet(out_parquet)
    return out_parquet

