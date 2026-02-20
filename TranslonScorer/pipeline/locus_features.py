from __future__ import annotations

from typing import Dict, List, Tuple

import polars as pl
from ..utils.logging import log_info


def _read_gtf(gtf_path: str) -> pl.DataFrame:
    """Lightweight GTF parser for features we need.

    Returns a Polars DataFrame with columns:
    chr, feature, start, end, strand, gene_id, transcript_id
    """
    df = pl.read_csv(
        gtf_path,
        has_header=False,
        separator="\t",
        comment_char="#",
        columns=[
            "column_1",  # seqname
            "column_3",  # feature
            "column_4",  # start
            "column_5",  # end
            "column_7",  # strand
            "column_9",  # attributes
        ],
        dtypes={
            "column_1": pl.Utf8,
            "column_3": pl.Utf8,
            "column_4": pl.Int64,
            "column_5": pl.Int64,
            "column_7": pl.Utf8,
            "column_9": pl.Utf8,
        },
        ignore_errors=True,
    ).rename(
        {
            "column_1": "chr",
            "column_3": "feature",
            "column_4": "start",
            "column_5": "end",
            "column_7": "strand",
            "column_9": "attributes",
        }
    )

    # Normalize chromosome names: allow 'chr' or bare names; keep original text
    # Extract gene_id and transcript_id
    df = df.with_columns(
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"').alias("gene_id"),
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"').alias("transcript_id"),
    ).drop("attributes")
    return df


def build_locus_features(gtf_path: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Build locus (gene) feature table and transcript→feature mappings.

    - exon_chunk features are derived by splitting exonic intervals on every transcript exon boundary within a gene
    - junction features are annotated exon-adjacency pairs within transcripts
    - TIS/TTS features from CDS starts/stops

    Returns:
        features_df: rows for each feature with genomic coords and membership
        map_df: rows mapping each transcript to an ordered sequence of feature_ids with transcript ranges
    """
    gtf = _read_gtf(gtf_path)

    exons = gtf.filter(pl.col("feature") == "exon").select(
        ["chr", "start", "end", "strand", "gene_id", "transcript_id"]
    )
    cds = gtf.filter(pl.col("feature") == "CDS").select(
        ["chr", "start", "end", "strand", "gene_id", "transcript_id"]
    )

    # For each gene, collect all exon boundary breakpoints
    def _chunks_for_gene(df_gene: pl.DataFrame) -> List[tuple[int, int]]:
        bounds: List[int] = []
        for s, e in zip(df_gene["start"], df_gene["end"]):
            bounds.append(int(s))
            bounds.append(int(e))
        if not bounds:
            return []
        bounds = sorted(set(bounds))
        chunks: List[tuple[int, int]] = []
        for i in range(len(bounds) - 1):
            a, b = bounds[i], bounds[i + 1]
            if a < b:
                chunks.append((a, b))
        return chunks

    features: List[dict] = []
    tmap_rows: List[dict] = []

    for (gene_id, chr_, strand), df_gene in (
        exons.group_by(["gene_id", "chr", "strand"]).agg(
            [pl.col("start"), pl.col("end"), pl.col("transcript_id")]
        ).iter_rows(named=True)
    ):
        # Expand row lists back to a tidy exon table
        starts = df_gene["start"]
        ends = df_gene["end"]
        trans = df_gene["transcript_id"]
        exon_tbl = pl.DataFrame({"start": starts, "end": ends, "transcript_id": trans})
        # Build exon chunks
        chunks = _chunks_for_gene(pl.DataFrame({"start": starts, "end": ends}))
        # Determine which chunks are exonic in any transcript (overlap test)
        chunk_ids: List[str] = []
        for ci, (a, b) in enumerate(chunks):
            # Exonic if covered by any exon interval
            covered = exon_tbl.filter((pl.col("start") < b) & (pl.col("end") > a))
            if covered.height == 0:
                continue
            fid = f"{gene_id}|chunk|{a}-{b}"
            chunk_ids.append(fid)
            features.append(
                {
                    "feature_id": fid,
                    "feature_type": "exon_chunk",
                    "locus_id": gene_id,
                    "chr": chr_,
                    "start": int(a),
                    "end": int(b),
                    "strand": strand,
                }
            )

        # Junction features from transcript exon adjacency
        # Reconstruct exons per transcript for this gene
        per_tx = (
            exon_tbl.group_by("transcript_id")
            .agg([pl.col("start").sort(), pl.col("end").sort()])
            .iter_rows(named=True)
        )
        for row in per_tx:
            tx = row["transcript_id"]
            s_list = row["start"]
            e_list = row["end"]
            # Produce ordered mapping of tx across chunks and junctions
            # Map exonic segments to chunk ids in order
            path: List[str] = []
            tx_tr_starts: List[int] = []
            tx_tr_ends: List[int] = []
            tx_tr_pos: List[int | None] = []
            tr_offset = 0
            for s, e in zip(s_list, e_list):
                # Identify chunks fully within this exon
                for (a, b) in chunks:
                    if a >= s and b <= e:
                        fid = f"{gene_id}|chunk|{a}-{b}"
                        if fid in chunk_ids:
                            path.append(fid)
                            tx_tr_starts.append(tr_offset + (a - s))
                            tx_tr_ends.append(tr_offset + (b - s))
                            tx_tr_pos.append(None)
                tr_offset += int(e) - int(s)
            # Insert junction features between consecutive exons for this transcript
            for i in range(len(s_list) - 1):
                donor = int(s_list[i + 0] if strand == "+" else e_list[i])
                acceptor = int(e_list[i] if strand == "+" else s_list[i + 1])
                jf = f"{gene_id}|junc|{donor}-{acceptor}"
                features.append(
                    {
                        "feature_id": jf,
                        "feature_type": "junction",
                        "locus_id": gene_id,
                        "chr": chr_,
                        "donor_pos": donor,
                        "acceptor_pos": acceptor,
                        "strand": strand,
                    }
                )
                path.append(jf)
                # Junction has zero length in transcript space; we store the position index as end of upstream chunk
                tx_tr_starts.append(None)
                tx_tr_ends.append(None)
                tx_tr_pos.append(None)

            # TIS/TTS from CDS
            cds_tx = cds.filter(pl.col("transcript_id") == tx)
            if cds_tx.height > 0:
                tis_pos = int(cds_tx["start"].min()) if strand == "+" else int(cds_tx["end"].max())
                tts_pos = int(cds_tx["end"].max()) if strand == "+" else int(cds_tx["start"].min())
                tis_id = f"{gene_id}|TIS|{tis_pos}"
                tts_id = f"{gene_id}|TTS|{tts_pos}"
                features.extend(
                    [
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
                    ]
                )
                # Compute transcript coordinate for TIS/TTS by walking exons
                def genomic_to_tran_pos(gpos: int) -> int | None:
                    tpos = 0
                    for s, e in zip(s_list, e_list):
                        s_i, e_i = int(s), int(e)
                        if s_i <= gpos < e_i:
                            return tpos + (gpos - s_i)
                        tpos += (e_i - s_i)
                    return None
                tis_tr = genomic_to_tran_pos(tis_pos)
                tts_tr = genomic_to_tran_pos(tts_pos)
                # Prepend/append TIS/TTS with transcript positions
                path = [tis_id] + path + [tts_id]
                tx_tr_starts = [None] + tx_tr_starts + [None]
                tx_tr_ends = [None] + tx_tr_ends + [None]
                tx_tr_pos = [tis_tr] + tx_tr_pos + [tts_tr]

            tmap_rows.append(
                {
                    "locus_id": gene_id,
                    "transcript_id": tx,
                    "feature_chain": path,
                    "tran_ranges_start": tx_tr_starts,
                    "tran_ranges_end": tx_tr_ends,
                    "tran_pos": tx_tr_pos,
                }
            )

    features_df = pl.from_dicts(features).unique(subset=["feature_id"], maintain_order=True)
    map_df = pl.from_dicts(tmap_rows)
    log_info(f"Built features: {features_df.height} rows; mappings: {map_df.height} rows")
    return features_df, map_df
