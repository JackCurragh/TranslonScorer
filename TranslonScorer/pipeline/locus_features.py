from __future__ import annotations

from typing import List, Tuple
import numpy as np

import pandas as pd
import polars as pl
import pyranges as pr
import click
import time
from ..utils.logging import log_info


def _read_gtf(gtf_path: str) -> pl.DataFrame:
    """Read GTF via PyRanges and return a Polars DataFrame
    with columns: chr, feature, start, end, strand, gene_id, transcript_id.
    """
    gr = pr.read_gtf(gtf_path)
    df = gr.df  # pandas DataFrame
    # Standard columns
    chrom = df["Chromosome"]
    start = df["Start"]
    end = df["End"]
    strand = df.get("Strand")
    feature = df.get("Feature")
    # Attributes (fallbacks if GTF uses different names)
    gene_col = next((c for c in ("gene_id", "gene", "gene_name", "Gene") if c in df.columns), None)
    tx_col = next((c for c in ("transcript_id", "transcript", "transcript_name", "Transcript") if c in df.columns), None)
    gene = df[gene_col] if gene_col else pd.Series([pd.NA] * len(df))
    tx = df[tx_col] if tx_col else pd.Series([pd.NA] * len(df))
    pdf = pd.DataFrame(
        {
            "chr": chrom,
            "feature": feature,
            "start": start,
            "end": end,
            "strand": strand,
            "gene_id": gene,
            "transcript_id": tx,
        }
    )
    return pl.from_pandas(pdf)


def build_locus_features(gtf_path: str, progress: bool = True) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Build locus (gene) feature table and transcript→feature mappings with PyRanges.

    - exon_chunk: maximal contiguous exonic intervals split at any exon boundary across transcripts (per gene/strand)
    - junction: exon adjacency within each transcript
    - TIS/TTS: from CDS blocks per transcript

    Returns
    -------
    features_df: Polars DataFrame with feature rows
    map_df: Polars DataFrame mapping transcript to ordered feature_chain with transcript-space ranges
    """
    t0=time.time()
    log_info("Reading annotation via PyRanges …")
    gtf = _read_gtf(gtf_path)
    log_info(f"Read annotation: {gtf.height:,} rows in {time.time()-t0:.2f}s")

    exons = gtf.filter(pl.col("feature") == "exon").select(
        ["chr", "start", "end", "strand", "gene_id", "transcript_id"]
    )
    cds = gtf.filter(pl.col("feature") == "CDS").select(
        ["chr", "start", "end", "strand", "gene_id", "transcript_id"]
    )

    # Utility: build chunk boundaries from all exon starts/ends within gene
    def _chunks_for_gene(df_gene: pl.DataFrame) -> List[Tuple[int, int]]:
        """Return maximal covered genomic spans for this gene by sweeping exon start/end events.

        This avoids per-chunk dataframe filtering by computing coverage once.
        """
        if df_gene.is_empty():
            return []
        starts = [int(x) for x in df_gene["start"].to_list()]
        ends = [int(x) for x in df_gene["end"].to_list()]
        if not starts:
            return []
        # Build event map
        ev = {}
        for s in starts:
            ev[s] = ev.get(s, 0) + 1
        for e in ends:
            ev[e] = ev.get(e, 0) - 1
        keys = sorted(ev.keys())
        spans: List[Tuple[int, int]] = []
        active = 0
        for i in range(len(keys) - 1):
            pos = keys[i]
            active += ev[pos]
            nxt = keys[i + 1]
            if active > 0 and pos < nxt:
                spans.append((pos, nxt))
        return spans

    features: List[dict] = []
    tmap_rows: List[dict] = []

    t1=time.time()
    # Precompute per-transcript exon arrays (sorted later by strand within each gene)
    tx_exons = (
        exons.group_by(["transcript_id", "strand"]).agg([
            pl.col("start").implode().alias("tx_starts"),
            pl.col("end").implode().alias("tx_ends"),
            pl.col("chr").first().alias("chr"),
            pl.col("gene_id").first().alias("gene_id"),
        ])
    )
    # Precompute TIS/TTS per transcript from CDS once
    cds_tx = (
        cds.group_by(["transcript_id", "strand"]).agg([
            pl.col("start").min().alias("cds_start_min"),
            pl.col("end").max().alias("cds_end_max"),
        ])
    )
    cds_lookup = { (r["transcript_id"], r["strand"]): (int(r["cds_start_min"]), int(r["cds_end_max"])) for r in cds_tx.iter_rows(named=True) }

    grouped = (
        tx_exons.group_by(["gene_id", "chr", "strand"]).agg([
            pl.col("transcript_id").implode().alias("transcripts"),
            pl.col("tx_starts").implode().alias("starts"),
            pl.col("tx_ends").implode().alias("ends"),
        ])
    )
    n_groups=grouped.height
    log_info(f"Index exons/CDS and building features for {n_groups:,} genes …")
    rows_iter = grouped.iter_rows(named=True)
    if progress:
        with click.progressbar(length=n_groups, label="Locus features") as bar:
            for row in rows_iter:
                bar.update(1)
                gene_id = row["gene_id"]
                chr_ = row["chr"]
                strand = row["strand"]
                starts = row["starts"]  # list[list[int]] per transcript
                ends = row["ends"]      # list[list[int]] per transcript
                trans = row["transcripts"]

                # Flatten all exons for chunk sweep
                all_s: List[int] = []
                all_e: List[int] = []
                for lst in starts:
                    all_s.extend(map(int, lst))
                for lst in ends:
                    all_e.extend(map(int, lst))
                exon_tbl = pl.DataFrame({"start": all_s, "end": all_e})

                # Exonic chunks from sweep-line
                chunks = _chunks_for_gene(exon_tbl.select(["start", "end"]))
                chunk_ids: List[str] = []
                for a, b in chunks:
                    fid = f"{gene_id}|chunk|{a}-{b}"
                    chunk_ids.append(fid)
                    features.append({
                        "feature_id": fid,
                        "feature_type": "exon_chunk",
                        "locus_id": gene_id,
                        "chr": chr_,
                        "start": int(a),
                        "end": int(b),
                        "strand": strand,
                    })

                # Transcript paths
                for tx, s_list, e_list in zip(trans, starts, ends):
                    s_list = list(map(int, s_list))
                    e_list = list(map(int, e_list))
                    # Order exons in transcript 5→3
                    if strand == "+":
                        order = sorted(range(len(s_list)), key=lambda i: s_list[i])
                    else:
                        order = sorted(range(len(s_list)), key=lambda i: s_list[i], reverse=True)
                    s_list = [s_list[i] for i in order]
                    e_list = [e_list[i] for i in order]

                    path: List[str] = []
                    tr_starts: List[int | None] = []
                    tr_ends: List[int | None] = []
                    tr_pos: List[int | None] = []
                    tr_offset = 0

                    for s, e in zip(s_list, e_list):
                        for (a, b) in chunks:
                            if a >= s and b <= e:
                                fid = f"{gene_id}|chunk|{a}-{b}"
                                if fid in chunk_ids:
                                    path.append(fid)
                                    if strand == "+":
                                        tr_starts.append(tr_offset + (a - s))
                                        tr_ends.append(tr_offset + (b - s))
                                    else:
                                        tr_starts.append(tr_offset + (e - b))
                                        tr_ends.append(tr_offset + (e - a))
                                    tr_pos.append(None)
                        tr_offset += (e - s)

                    # Junctions between consecutive exons
                    for i in range(len(s_list) - 1):
                        if strand == "+":
                            donor, acceptor = int(e_list[i]), int(e_list[i + 1])
                        else:
                            donor, acceptor = int(s_list[i]), int(s_list[i + 1])
                        jf = f"{gene_id}|junc|{donor}-{acceptor}"
                        features.append({
                            "feature_id": jf,
                            "feature_type": "junction",
                            "locus_id": gene_id,
                            "chr": chr_,
                            "donor_pos": donor,
                            "acceptor_pos": acceptor,
                            "strand": strand,
                        })
                        path.append(jf)
                        tr_starts.append(None)
                        tr_ends.append(None)
                        tr_pos.append(None)

                    # TIS/TTS from precomputed CDS
                    cds_pair = cds_lookup.get((tx, strand))
                    if cds_pair is not None:
                        if strand == "+":
                            tis_pos = int(cds_pair[0])
                            tts_pos = int(cds_pair[1])
                        else:
                            tis_pos = int(cds_pair[1])
                            tts_pos = int(cds_pair[0])
                        tis_id = f"{gene_id}|TIS|{tis_pos}"
                        tts_id = f"{gene_id}|TTS|{tts_pos}"
                        features.extend([
                            {
                                "feature_id": tis_id,
                                "feature_type": "TIS",
                                "locus_id": gene_id,
                                "chr": chr_,
                                "pos": tis_pos,
                                "strand": strand,
                            },
                            {
                                "feature_id": tts_id,
                                "feature_type": "TTS",
                                "locus_id": gene_id,
                                "chr": chr_,
                                "pos": tts_pos,
                                "strand": strand,
                            },
                        ])

                        def genomic_to_tran_pos(gpos: int) -> int | None:
                            tpos = 0
                            for s, e in zip(s_list, e_list):
                                if s <= gpos < e:
                                    if strand == "+":
                                        return tpos + (gpos - s)
                                    else:
                                        return tpos + (e - gpos - 1)
                                tpos += (e - s)
                            return None

                        tis_tr = genomic_to_tran_pos(tis_pos)
                        tts_tr = genomic_to_tran_pos(tts_pos)
                        path = [tis_id] + path + [tts_id]
                        tr_starts = [None] + tr_starts + [None]
                        tr_ends = [None] + tr_ends + [None]
                        tr_pos = [tis_tr] + tr_pos + [tts_tr]

                    tmap_rows.append({
                        "locus_id": gene_id,
                        "transcript_id": tx,
                        "feature_chain": path,
                        "tran_ranges_start": tr_starts,
                        "tran_ranges_end": tr_ends,
                        "tran_pos": tr_pos,
                    })
    else:
        for row in rows_iter:
            gene_id = row["gene_id"]
            chr_ = row["chr"]
            strand = row["strand"]
            starts = row["starts"]
            ends = row["ends"]
            trans = row["transcripts"]
            exon_tbl = pl.DataFrame({"start": starts, "end": ends, "transcript_id": trans})
            chunks = _chunks_for_gene(exon_tbl.select(["start", "end"]))
            chunk_ids: List[str] = []
            for a, b in chunks:
                covered = exon_tbl.filter((pl.col("start") < b) & (pl.col("end") > a))
                if covered.height == 0:
                    continue
                fid = f"{gene_id}|chunk|{a}-{b}"
                chunk_ids.append(fid)
                features.append({
                    "feature_id": fid,
                    "feature_type": "exon_chunk",
                    "locus_id": gene_id,
                    "chr": chr_,
                    "start": int(a),
                    "end": int(b),
                    "strand": strand,
                })
            per_tx = (
                exon_tbl.group_by("transcript_id")
                .agg([pl.col("start").sort(), pl.col("end").sort()])
                .iter_rows(named=True)
            )
            for rtx in per_tx:
                tx = rtx["transcript_id"]
                s_list = list(map(int, rtx["start"]))
                e_list = list(map(int, rtx["end"]))
                if strand == "+":
                    order = sorted(range(len(s_list)), key=lambda i: s_list[i])
                else:
                    order = sorted(range(len(s_list)), key=lambda i: s_list[i], reverse=True)
                s_list = [s_list[i] for i in order]
                e_list = [e_list[i] for i in order]
                path: List[str] = []
                tr_starts: List[int | None] = []
                tr_ends: List[int | None] = []
                tr_pos: List[int | None] = []
                tr_offset = 0
                for s, e in zip(s_list, e_list):
                    for (a, b) in chunks:
                        if a >= s and b <= e:
                            fid = f"{gene_id}|chunk|{a}-{b}"
                            if fid in chunk_ids:
                                path.append(fid)
                                if strand == "+":
                                    tr_starts.append(tr_offset + (a - s))
                                    tr_ends.append(tr_offset + (b - s))
                                else:
                                    tr_starts.append(tr_offset + (e - b))
                                    tr_ends.append(tr_offset + (e - a))
                                tr_pos.append(None)
                    tr_offset += (e - s)
                for i in range(len(s_list) - 1):
                    if strand == "+":
                        donor, acceptor = int(e_list[i]), int(s_list[i + 1])
                    else:
                        donor, acceptor = int(s_list[i]), int(e_list[i + 1])
                    jf = f"{gene_id}|junc|{donor}-{acceptor}"
                    features.append({
                        "feature_id": jf,
                        "feature_type": "junction",
                        "locus_id": gene_id,
                        "chr": chr_,
                        "donor_pos": donor,
                        "acceptor_pos": acceptor,
                        "strand": strand,
                    })
                    path.append(jf)
                    tr_starts.append(None)
                    tr_ends.append(None)
                    tr_pos.append(None)
                cds_tx = cds.filter(pl.col("transcript_id") == tx)
                if cds_tx.height > 0:
                    if strand == "+":
                        tis_pos = int(cds_tx["start"].min())
                        tts_pos = int(cds_tx["end"].max())
                    else:
                        tis_pos = int(cds_tx["end"].max())
                        tts_pos = int(cds_tx["start"].min())
                    tis_id = f"{gene_id}|TIS|{tis_pos}"
                    tts_id = f"{gene_id}|TTS|{tts_pos}"
                    features.extend([
                        {
                            "feature_id": tis_id,
                            "feature_type": "TIS",
                            "locus_id": gene_id,
                            "chr": chr_,
                            "pos": tis_pos,
                            "strand": strand,
                        },
                        {
                            "feature_id": tts_id,
                            "feature_type": "TTS",
                            "locus_id": gene_id,
                            "chr": chr_,
                            "pos": tts_pos,
                            "strand": strand,
                        },
                    ])
                    def genomic_to_tran_pos(gpos: int) -> int | None:
                        tpos = 0
                        for s, e in zip(s_list, e_list):
                            if s <= gpos < e:
                                if strand == "+":
                                    return tpos + (gpos - s)
                                else:
                                    return tpos + (e - gpos - 1)
                            tpos += (e - s)
                        return None
                    tis_tr = genomic_to_tran_pos(tis_pos)
                    tts_tr = genomic_to_tran_pos(tts_pos)
                    path = [tis_id] + path + [tts_id]
                    tr_starts = [None] + tr_starts + [None]
                    tr_ends = [None] + tr_ends + [None]
                    tr_pos = [tis_tr] + tr_pos + [tts_tr]
                tmap_rows.append({
                    "locus_id": gene_id,
                    "transcript_id": tx,
                    "feature_chain": path,
                    "tran_ranges_start": tr_starts,
                    "tran_ranges_end": tr_ends,
                    "tran_pos": tr_pos,
                })
        gene_id = row["gene_id"]
        chr_ = row["chr"]
        strand = row["strand"]
        starts = row["starts"]
        ends = row["ends"]
        trans = row["transcripts"]
        # Flatten all exon starts/ends for chunk sweep
        all_s = []
        all_e = []
        for lst in starts:
            all_s.extend(map(int, lst))
        for lst in ends:
            all_e.extend(map(int, lst))
        exon_tbl = pl.DataFrame({"start": all_s, "end": all_e})

        # Exonic chunks (split at any exon boundary, keep spans covered by at least one exon)
        chunks = _chunks_for_gene(exon_tbl.select(["start", "end"]))
        chunk_ids: List[str] = []
        for a, b in chunks:
            fid = f"{gene_id}|chunk|{a}-{b}"
            chunk_ids.append(fid)
            features.append({
                "feature_id": fid,
                "feature_type": "exon_chunk",
                "locus_id": gene_id,
                "chr": chr_,
                "start": int(a),
                "end": int(b),
                "strand": strand,
            })

        # Per-transcript path: exonic chunks in transcript order + junctions + TIS/TTS
        # Iterate transcripts listed in this gene
        for tx, s_list, e_list in zip(trans, starts, ends):
            s_list = list(map(int, s_list))
            e_list = list(map(int, e_list))
            # Order exons in transcript 5→3
            if strand == "+":
                order = sorted(range(len(s_list)), key=lambda i: s_list[i])
            else:
                order = sorted(range(len(s_list)), key=lambda i: s_list[i], reverse=True)
            s_list = [s_list[i] for i in order]
            e_list = [e_list[i] for i in order]

            path: List[str] = []
            tr_starts: List[int | None] = []
            tr_ends: List[int | None] = []
            tr_pos: List[int | None] = []
            tr_offset = 0

            for s, e in zip(s_list, e_list):
                for (a, b) in chunks:
                    if a >= s and b <= e:
                        fid = f"{gene_id}|chunk|{a}-{b}"
                        if fid in chunk_ids:
                            path.append(fid)
                            if strand == "+":
                                tr_starts.append(tr_offset + (a - s))
                                tr_ends.append(tr_offset + (b - s))
                            else:
                                # minus strand: transcript coords run from exon end to start
                                tr_starts.append(tr_offset + (e - b))
                                tr_ends.append(tr_offset + (e - a))
                            tr_pos.append(None)
                tr_offset += (e - s)

            # Junctions between consecutive exons
            for i in range(len(s_list) - 1):
                if strand == "+":
                    donor, acceptor = int(e_list[i]), int(s_list[i + 1])
                else:
                    donor, acceptor = int(s_list[i]), int(e_list[i + 1])
                jf = f"{gene_id}|junc|{donor}-{acceptor}"
                features.append({
                    "feature_id": jf,
                    "feature_type": "junction",
                    "locus_id": gene_id,
                    "chr": chr_,
                    "donor_pos": donor,
                    "acceptor_pos": acceptor,
                    "strand": strand,
                })
                path.append(jf)
                tr_starts.append(None)
                tr_ends.append(None)
                tr_pos.append(None)

            # TIS/TTS from CDS (if available for this transcript)
            cds_pair = cds_lookup.get((tx, strand))
            if cds_pair is not None:
                if strand == "+":
                    tis_pos = int(cds_pair[0])
                    tts_pos = int(cds_pair[1])
                else:
                    tis_pos = int(cds_pair[1])
                    tts_pos = int(cds_pair[0])
                tis_id = f"{gene_id}|TIS|{tis_pos}"
                tts_id = f"{gene_id}|TTS|{tts_pos}"
                features.extend([
                    {
                        "feature_id": tis_id,
                        "feature_type": "TIS",
                        "locus_id": gene_id,
                        "chr": chr_,
                        "pos": tis_pos,
                        "strand": strand,
                    },
                    {
                        "feature_id": tts_id,
                        "feature_type": "TTS",
                        "locus_id": gene_id,
                        "chr": chr_,
                        "pos": tts_pos,
                        "strand": strand,
                    },
                ])

                def genomic_to_tran_pos(gpos: int) -> int | None:
                    tpos = 0
                    for s, e in zip(s_list, e_list):
                        if s <= gpos < e:
                            if strand == "+":
                                return tpos + (gpos - s)
                            else:
                                return tpos + (e - gpos - 1)
                        tpos += (e - s)
                    return None

                tis_tr = genomic_to_tran_pos(tis_pos)
                tts_tr = genomic_to_tran_pos(tts_pos)
                path = [tis_id] + path + [tts_id]
                tr_starts = [None] + tr_starts + [None]
                tr_ends = [None] + tr_ends + [None]
                tr_pos = [tis_tr] + tr_pos + [tts_tr]

            tmap_rows.append({
                "locus_id": gene_id,
                "transcript_id": tx,
                "feature_chain": path,
                "tran_ranges_start": tr_starts,
                "tran_ranges_end": tr_ends,
                "tran_pos": tr_pos,
            })

    features_df = pl.from_dicts(features).unique(subset=["feature_id"], maintain_order=True)
    map_df = pl.from_dicts(tmap_rows)
    log_info(f"Built features: {features_df.height} rows; mappings: {map_df.height} rows in {time.time()-t1:.2f}s (total {time.time()-t0:.2f}s)")
    return features_df, map_df
