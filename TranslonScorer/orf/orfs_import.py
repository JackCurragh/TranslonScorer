from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Iterable, Set

import pandas as pd
import polars as pl
from ..utils.logging import log_info
import click
import re
import os

@dataclass
class TxModel:
    chrom: str
    strand: str
    exons: List[Tuple[int, int]]  # 5'->3' transcript order, half-open [start,end)
    cumlens: List[int]            # cumulative transcript length at exon starts
    length: int
    gene_id: Optional[str]

def _expand_chrom_filter(chroms: Set[str]) -> Set[str]:
    out: Set[str] = set()
    for c in chroms:
        s = str(c)
        out.add(s)
        if s.startswith('chr'):
            core = s[3:]
            out.add(core)
            if core == 'M':
                out.add('MT')
        else:
            out.add('chr' + s)
            if s == 'MT':
                out.add('chrM')
            if s == 'M':
                out.add('MT'); out.add('chrM')
    return out

def _norm_chr(s: str) -> str:
    ss = str(s)
    if ss.startswith('chr'):
        return ss
    if ss in ('MT', 'M', 'Mt', 'mt'):
        return 'chrM'
    return 'chr' + ss


def _scan_gtf_exons(gtf_path: str, chrom_filter: Optional[Set[str]] = None) -> Iterable[Dict[str, object]]:
    """Stream exon rows from a GTF using pandas chunks (low-memory, comment-aware).

    Yields dicts with keys: Chromosome, Start, End, Strand, transcript_id, gene_id.
    """
    cols = [
        "Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes",
    ]
    usecols = [0, 2, 3, 4, 6, 8]
    log_info("Scanning GTF exons (chunked pandas)…")
    tx_re = re.compile(r'(?:^|;)\s*transcript_id\s+"?([^";]+)"?')
    gene_re = re.compile(r'(?:^|;)\s*gene_id\s+"?([^";]+)"?')
    for chunk in pd.read_csv(
        gtf_path,
        sep="\t",
        header=None,
        names=cols,
        usecols=["Chromosome", "Feature", "Start", "End", "Strand", "Attributes"],
        comment="#",
        dtype={
            "Chromosome": str,
            "Feature": str,
            "Start": int,
            "End": int,
            "Strand": str,
            "Attributes": str,
        },
        chunksize=200_000,
        engine="c",
        on_bad_lines="skip",
    ):
        ex = chunk[chunk["Feature"] == "exon"].copy()
        if chrom_filter is not None:
            ex = ex[ex["Chromosome"].astype(str).isin(_expand_chrom_filter(chrom_filter))]
        if ex.empty:
            continue
        # Extract transcript_id/gene_id
        ex["transcript_id"] = ex["Attributes"].str.extract(tx_re)
        ex["gene_id"] = ex["Attributes"].str.extract(gene_re)
        ex = ex.drop(columns=["Attributes", "Feature", "Source", "Score", "Frame"], errors="ignore")
        for row in ex.itertuples(index=False):
            yield {
                "Chromosome": row.Chromosome,
                "Start": int(row.Start),
                "End": int(row.End),
                "Strand": row.Strand,
                "transcript_id": row.transcript_id,
                "gene_id": row.gene_id,
            }


def _build_tx_models_streaming(gtf_path: str, chroms: Set[str]) -> Dict[str, TxModel]:
    tx_models: Dict[str, TxModel] = {}
    exons_by_tx: Dict[str, List[Tuple[int, int]]] = {}
    meta_by_tx: Dict[str, Tuple[str, str, Optional[str]]] = {}
    for r in _scan_gtf_exons(gtf_path, chrom_filter=chroms):
        tx = r.get("transcript_id")
        if tx is None:
            # Skip entries lacking transcript_id (non-standard GTF)
            continue
        chrom = str(r["Chromosome"]) if r.get("Chromosome") is not None else None
        strand = str(r.get("Strand")) if r.get("Strand") is not None else None
        s = int(r["Start"]) - 1  # GTF is 1-based inclusive; convert to 0-based half-open
        e = int(r["End"])       # end becomes exclusive
        if chrom is None or strand is None:
            continue
        exons_by_tx.setdefault(str(tx), []).append((s, e))
        if str(tx) not in meta_by_tx:
            meta_by_tx[str(tx)] = (_norm_chr(chrom), strand, r.get("gene_id"))

    for tx, exons in exons_by_tx.items():
        chrom, strand, gene_id = meta_by_tx.get(tx, (None, None, None))
        if chrom is None or strand is None:
            continue
        # Sort to transcript 5'→3' order
        if strand == '+':
            exons = sorted(exons, key=lambda t: t[0])
        else:
            exons = sorted(exons, key=lambda t: t[0], reverse=True)
        cum, acc = [], 0
        for s, e in exons:
            cum.append(acc)
            acc += (e - s)
        tx_models[tx] = TxModel(chrom, strand, exons, cum, acc, str(gene_id) if gene_id is not None else None)
    return tx_models

def _bed12_to_blocks(row: pd.Series) -> List[Tuple[int, int]]:
    cs = int(row[1]); ce = int(row[2])
    bc = int(row[9]) if not pd.isna(row[9]) else 1
    sizes = [int(x) for x in str(row[10]).strip(',').split(',') if x]
    starts = [int(x) for x in str(row[11]).strip(',').split(',') if x]
    if bc != len(sizes) or bc != len(starts) or bc == 0:
        return [(cs, ce)]
    out = []
    for i in range(bc):
        s = cs + starts[i]
        e = s + sizes[i]
        if e > s:
            out.append((s, e))
    return out

def _junctions_from_blocks(blocks: List[Tuple[int,int]], strand: str) -> List[Tuple[int,int]]:
    if len(blocks) <= 1:
        return []
    juncs = []
    if strand == '+':
        for i in range(len(blocks)-1):
            juncs.append((blocks[i][1], blocks[i+1][0]))
    else:
        for i in range(len(blocks)-1):
            juncs.append((blocks[i][0], blocks[i+1][1]))
    return juncs

def _overlap_len(a: Tuple[int,int], b: Tuple[int,int]) -> int:
    s = max(a[0], b[0]); e = min(a[1], b[1])
    return max(0, e - s)

def _map_genomic_to_tran(tx: TxModel, gpos: int) -> Optional[int]:
    for idx, (s, e) in enumerate(tx.exons):
        if s <= gpos < e:
            if tx.strand == '+':
                return tx.cumlens[idx] + (gpos - s)
            else:
                return tx.cumlens[idx] + (e - gpos - 1)
    return None

def import_bed12(
    bed12_path: str,
    gtf_path: str,
    out_parquet: str,
    assign_policy: str = 'best',
    require_junction_match: bool = True,
    progress: bool = True,
) -> pl.DataFrame:
    """Import ORFs from BED12, map to transcripts, and write canonical ORFs (Parquet)."""
    # Read BED12 first to determine chromosomes to load from GTF (reduces memory)
    log_info("Reading BED12 candidates…")
    bed = pd.read_csv(
        bed12_path, sep='\t', header=None, comment='#', dtype={0:str,1:int,2:int,3:str,5:str}, engine='python'
    ).dropna(subset=[0,1,2,5])
    log_info(f"Loaded BED12: {len(bed):,} rows")
    chroms_needed: Set[str] = set(map(str, bed[0].unique().tolist()))
    log_info(f"Building transcript models for {len(chroms_needed)} chromosome(s)…")
    tx_models = _build_tx_models_streaming(gtf_path, chroms_needed)
    log_info(f"Built models for {len(tx_models):,} transcripts")
    by_chr_strand: Dict[Tuple[str,str], List[str]] = {}
    tx_bounds: Dict[str, Tuple[int,int]] = {}
    for tid, tx in tx_models.items():
        by_chr_strand.setdefault((tx.chrom, tx.strand), []).append(tid)
        tx_bounds[tid] = (min(s for s,_ in tx.exons), max(e for _,e in tx.exons))
    out_rows = []
    log_info("Mapping ORFs to transcripts…")
    row_iter = bed.iterrows()
    if progress:
        with click.progressbar(length=len(bed), label="Map ORFs") as bar:
            for idx, row in row_iter:
                bar.update(1)
                chrom = _norm_chr(str(row[0])); start = int(row[1]); end = int(row[2])
                name = str(row[3]) if not pd.isna(row[3]) else f'orf_{idx}'
                strand = str(row[5])
                blocks = _bed12_to_blocks(row)
                if not blocks:
                    continue
                juncs = _junctions_from_blocks(blocks, strand)
                length_nt = sum(b2-b1 for b1,b2 in blocks)
                candidates = by_chr_strand.get((chrom, strand), [])
                scored: List[Tuple[int,int,str]] = []
                for tid in candidates:
                    tx = tx_models[tid]
                    tb = tx_bounds[tid]
                    if _overlap_len((start,end), tb) == 0:
                        continue
                    ov = 0
                    for b in blocks:
                        for ex in tx.exons:
                            ov += _overlap_len(b, ex)
                    tj = _junctions_from_blocks(tx.exons, strand)
                    jm = sum(1 for jj in juncs if jj in tj)
                    if require_junction_match and len(juncs) > 0 and jm != len(juncs):
                        continue
                    if ov == 0:
                        continue
                    scored.append((jm, ov, tid))
                if not scored:
                    out_rows.append({
                        'orf_id': f'bed12:{name}:{chrom}:{start}-{end}:{strand}',
                        'source': 'bed12', 'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                        'block_starts': [b[0] for b in blocks], 'block_ends': [b[1] for b in blocks],
                        'length_nt': length_nt, 'tran_id': None, 'start_pos_tran': None, 'stop_pos_tran': None,
                        'transcript_ambiguity': True,
                    })
                    continue
                if assign_policy == 'best':
                    scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
                    keep = [scored[0]]
                elif assign_policy == 'all':
                    keep = scored
                else:
                    keep = scored
                for jm, ov, tid in keep:
                    tx = tx_models[tid]
                    gstart = blocks[0][0]
                    gend_last = blocks[-1][1] - 1
                    t_start = _map_genomic_to_tran(tx, gstart)
                    if t_start is not None:
                        t_end = t_start + length_nt
                    else:
                        te_last = _map_genomic_to_tran(tx, gend_last)
                        t_end = te_last + 1 if te_last is not None else None
                    out_rows.append({
                        'orf_id': f'bed12:{name}:{chrom}:{start}-{end}:{strand}:{tid}',
                        'source': 'bed12', 'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                        'block_starts': [b[0] for b in blocks], 'block_ends': [b[1] for b in blocks],
                        'length_nt': length_nt, 'tran_id': tid,
                        'start_pos_tran': t_start, 'stop_pos_tran': t_end,
                        'transcript_ambiguity': (assign_policy=='all' and len(keep)>1),
                    })
    else:
        for idx, row in row_iter:
            chrom = _norm_chr(str(row[0])); start = int(row[1]); end = int(row[2])
            name = str(row[3]) if not pd.isna(row[3]) else f'orf_{idx}'
            strand = str(row[5])
            blocks = _bed12_to_blocks(row)
            if not blocks:
                continue
            juncs = _junctions_from_blocks(blocks, strand)
            length_nt = sum(b2-b1 for b1,b2 in blocks)
            candidates = by_chr_strand.get((chrom, strand), [])
            scored: List[Tuple[int,int,str]] = []
            for tid in candidates:
                tx = tx_models[tid]
                tb = tx_bounds[tid]
                if _overlap_len((start,end), tb) == 0:
                    continue
                ov = 0
                for b in blocks:
                    for ex in tx.exons:
                        ov += _overlap_len(b, ex)
                tj = _junctions_from_blocks(tx.exons, strand)
                jm = sum(1 for jj in juncs if jj in tj)
                if require_junction_match and len(juncs) > 0 and jm != len(juncs):
                    continue
                if ov == 0:
                    continue
                scored.append((jm, ov, tid))
            if not scored:
                out_rows.append({
                    'orf_id': f'bed12:{name}:{chrom}:{start}-{end}:{strand}',
                    'source': 'bed12', 'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                    'block_starts': [b[0] for b in blocks], 'block_ends': [b[1] for b in blocks],
                    'length_nt': length_nt, 'tran_id': None, 'start_pos_tran': None, 'stop_pos_tran': None,
                    'transcript_ambiguity': True,
                })
                continue
            if assign_policy == 'best':
                scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
                keep = [scored[0]]
            elif assign_policy == 'all':
                keep = scored
            else:
                keep = scored
            for jm, ov, tid in keep:
                tx = tx_models[tid]
                gstart = blocks[0][0]
                gend_last = blocks[-1][1] - 1
                t_start = _map_genomic_to_tran(tx, gstart)
                if t_start is not None:
                    t_end = t_start + length_nt
                else:
                    te_last = _map_genomic_to_tran(tx, gend_last)
                    t_end = te_last + 1 if te_last is not None else None
                out_rows.append({
                    'orf_id': f'bed12:{name}:{chrom}:{start}-{end}:{strand}:{tid}',
                    'source': 'bed12', 'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                    'block_starts': [b[0] for b in blocks], 'block_ends': [b[1] for b in blocks],
                    'length_nt': length_nt, 'tran_id': tid,
                    'start_pos_tran': t_start, 'stop_pos_tran': t_end,
                    'transcript_ambiguity': (assign_policy=='all' and len(keep)>1),
                })
    df = pl.from_pandas(pd.DataFrame(out_rows))
    os.makedirs(os.path.dirname(out_parquet) or '.', exist_ok=True)
    df.write_parquet(out_parquet)
    return df
