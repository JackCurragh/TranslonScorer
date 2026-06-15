from __future__ import annotations

from typing import List, Dict, Any, Tuple, Optional
import click

import polars as pl


def _slice_orf_chain(
    path: List[str],
    tr_starts: List[int | None],
    tr_ends: List[int | None],
    tr_pos: List[int | None],
    orf_start: int,
    orf_end: int,
    types_by_id: Dict[str, str] | None = None,
):
    """Return per-ORF feature chain and slices relative to ORF.

    - For range features (exon_chunk): slice_start/slice_end are ORF-relative [start, end).
    - For point features (TIS/TTS): event_pos is ORF-relative.
    - Junctions are included if both adjacent chunks overlap the ORF; they carry no slice range.
    """
    n = len(path)
    out_ids: List[str] = []
    out_types: List[str] = []
    slice_start: List[int | None] = []
    slice_end: List[int | None] = []
    event_pos: List[int | None] = []

    # Mark chunk indices that overlap the ORF
    chunk_overlap: Dict[int, bool] = {}
    for i in range(n):
        ts, te = tr_starts[i], tr_ends[i]
        if ts is None or te is None:
            continue
        if te <= orf_start or ts >= orf_end:
            continue
        chunk_overlap[i] = True

    def _append(idx: int, s: int | None, e: int | None, p: int | None):
        out_ids.append(path[idx])
        ftype = types_by_id.get(path[idx], "unknown") if types_by_id else "unknown"
        out_types.append(ftype)
        slice_start.append(s)
        slice_end.append(e)
        event_pos.append(p)

    i = 0
    while i < n:
        ts, te = tr_starts[i], tr_ends[i]
        pos = tr_pos[i]
        fid = path[i]
        # Range feature
        if ts is not None and te is not None:
            if not (te <= orf_start or ts >= orf_end):
                s = max(ts, orf_start) - orf_start
                e = min(te, orf_end) - orf_start
                _append(i, int(s), int(e), None)
        else:
            # Point feature (TIS/TTS) with pos recorded
            if pos is not None:
                if orf_start <= pos < orf_end:
                    _append(i, None, None, int(pos - orf_start))
            else:
                # Likely a junction (between chunks). Include only if both adjacent chunks overlap ORF.
                left_ok = chunk_overlap.get(i - 1, False)
                right_ok = chunk_overlap.get(i + 1, False)
                if left_ok and right_ok:
                    _append(i, None, None, None)
        i += 1

    return out_ids, out_types, slice_start, slice_end, event_pos


def _build_composites(
    fr: Dict[str, Any],
    orf_start: int,
    orf_end: int,
    types_by_id: Dict[str, str],
    tis_range: int = 15,
    flank_nt: int = 60,
) -> Tuple[List[str], List[str], List[Optional[int]], List[Optional[int]]]:
    """Construct composite parts around an ORF: upstream chunk, TIS window, ORF chunks + junctions, TTS window, downstream chunk.

    Returns arrays (comp_kind, comp_feature_id, comp_slice_start, comp_slice_end) with ORF-relative coordinates.
    """
    chain = fr["feature_chain"]
    trs = fr["tran_ranges_start"]
    tre = fr["tran_ranges_end"]
    tpos = fr.get("tran_pos") or [None] * len(chain)

    def find_chunk_at(pos: int, right_edge: bool = False) -> Optional[int]:
        for i, (a, b) in enumerate(zip(trs, tre)):
            if a is None or b is None:
                continue
            if not right_edge and (a <= pos < b):
                return i
            if right_edge and (a < pos <= b):
                return i
        return None

    comp_kind: List[str] = []
    comp_id: List[str] = []
    comp_s: List[Optional[int]] = []
    comp_e: List[Optional[int]] = []

    # Upstream chunk window
    si = find_chunk_at(orf_start, right_edge=False)
    if si is not None and types_by_id.get(chain[si]) == "exon_chunk":
        a, b = trs[si], tre[si]
        if a is not None and b is not None:
            us = max(a, orf_start - flank_nt)
            ue = orf_start
            if ue > us:
                comp_kind.append("upstream_chunk")
                comp_id.append(chain[si])
                comp_s.append(int(us - orf_start))
                comp_e.append(int(ue - orf_start))

    # TIS window around ORF start
    comp_kind.append("TIS_region")
    # Prefer actual TIS feature id if present in chain; fall back to literal
    tis_fid = next(
        (fid for fid, ty in zip(chain, (types_by_id.get(x, "") for x in chain)) if ty == "TIS"),
        "TIS",
    )
    comp_id.append(tis_fid)
    comp_s.append(-int(tis_range))
    comp_e.append(int(tis_range))

    # ORF chunks + interleaving junctions
    prev_in = False
    for i, fid in enumerate(chain):
        ftype = types_by_id.get(fid, "")
        if ftype == "exon_chunk":
            a, b = trs[i], tre[i]
            if a is None or b is None:
                prev_in = False
                continue
            if not (b <= orf_start or a >= orf_end):
                rs = max(a, orf_start) - orf_start
                re = min(b, orf_end) - orf_start
                comp_kind.append("orf_chunk")
                comp_id.append(fid)
                comp_s.append(int(rs))
                comp_e.append(int(re))
                prev_in = True
            else:
                prev_in = False
        elif ftype == "junction":
            if prev_in and (i + 1) < len(chain):
                a2, b2 = trs[i + 1], tre[i + 1]
                if a2 is not None and b2 is not None and not (b2 <= orf_start or a2 >= orf_end):
                    comp_kind.append("junction")
                    comp_id.append(fid)
                    comp_s.append(None)
                    comp_e.append(None)
            prev_in = False
        else:
            prev_in = False

    # TTS window around ORF stop
    comp_kind.append("TTS_region")
    tts_fid = next(
        (fid for fid, ty in zip(chain, (types_by_id.get(x, "") for x in chain)) if ty == "TTS"),
        "TTS",
    )
    comp_id.append(tts_fid)
    comp_s.append(int((orf_end - orf_start) - tis_range))
    comp_e.append(int((orf_end - orf_start) + tis_range))

    # Downstream chunk window
    ei = find_chunk_at(orf_end, right_edge=True)
    if ei is not None and types_by_id.get(chain[ei]) == "exon_chunk":
        a, b = trs[ei], tre[ei]
        if a is not None and b is not None:
            ds = orf_end
            de = min(b, orf_end + flank_nt)
            if de > ds:
                comp_kind.append("downstream_chunk")
                comp_id.append(chain[ei])
                comp_s.append(int(ds - orf_start))
                comp_e.append(int(de - orf_start))

    return comp_kind, comp_id, comp_s, comp_e


def map_orfs(
    orfs_parquet: str,
    feature_map_parquet: str,
    features_parquet: str,
    out_parquet: str,
    progress: bool = True,
    tis_range: int = 15,
    flank_nt: int = 60,
):
    """Map ORFs (transcript-space) to per-ORF feature chains and slices.

    Inputs
    ------
    - orfs_parquet: contains orf_id, tran_id, start_pos_tran, stop_pos_tran.
    - feature_map_parquet: per-transcript feature_chain with tran_ranges_* and tran_pos.
    - features_parquet: feature_id → feature_type (exon_chunk/junction/TIS/TTS) and locus metadata.
    """
    orfs = pl.read_parquet(orfs_parquet)
    fmap = pl.read_parquet(feature_map_parquet)
    feats = pl.read_parquet(features_parquet)

    # Build lookup: transcript_id -> row with arrays
    fmap_lookup: Dict[str, Any] = {}
    for r in fmap.iter_rows(named=True):
        fmap_lookup[str(r["transcript_id"])] = r

    # Feature type map
    types_by_id: Dict[str, str] = {
        str(fid): str(ft) for fid, ft in feats.select(["feature_id", "feature_type"]).iter_rows()
    }
    locus_by_id: Dict[str, str] = {
        str(fid): str(lid) for fid, lid in feats.select(["feature_id", "locus_id"]).iter_rows()
    }

    out_rows: List[Dict[str, Any]] = []
    rows_iter = orfs.iter_rows(named=True)
    if progress:
        with click.progressbar(length=orfs.height, label="Map ORFs→features") as bar:
            for r in rows_iter:
                bar.update(1)
                tid = r.get("tran_id") or r.get("transcript_id")
                s = r.get("start_pos_tran")
                e = r.get("stop_pos_tran")
                if tid is None or s is None or e is None:
                    continue
                s = int(s)
                e = int(e)
                if e <= s:
                    continue
                fr = fmap_lookup.get(str(tid))
                if fr is None:
                    continue
                ids, types, ss, ee, pos = _slice_orf_chain(
                    fr["feature_chain"],
                    fr["tran_ranges_start"],
                    fr["tran_ranges_end"],
                    fr["tran_pos"],
                    s,
                    e,
                    types_by_id,
                )
                if not ids:
                    continue
                locus = None
                for fid in ids:
                    locus = locus_by_id.get(fid)
                    if locus:
                        break
                # Composites
                ck, cid, cs, ce = _build_composites(
                    fr, s, e, types_by_id, tis_range=tis_range, flank_nt=flank_nt
                )
                out_rows.append(
                    {
                        "orf_id": r.get("orf_id"),
                        "tran_id": str(tid),
                        "start_pos_tran": s,
                        "stop_pos_tran": e,
                        "locus_id": locus,
                        "feature_chain": ids,
                        "feature_types": types,
                        "slice_start": ss,
                        "slice_end": ee,
                        "event_pos": pos,
                        "comp_kind": ck,
                        "comp_feature_id": cid,
                        "comp_slice_start": cs,
                        "comp_slice_end": ce,
                    }
                )
    else:
        for r in rows_iter:
            tid = r.get("tran_id") or r.get("transcript_id")
            s = r.get("start_pos_tran")
            e = r.get("stop_pos_tran")
            if tid is None or s is None or e is None:
                continue
            s = int(s)
            e = int(e)
            if e <= s:
                continue
            fr = fmap_lookup.get(str(tid))
            if fr is None:
                continue
            ids, types, ss, ee, pos = _slice_orf_chain(
                fr["feature_chain"],
                fr["tran_ranges_start"],
                fr["tran_ranges_end"],
                fr["tran_pos"],
                s,
                e,
                types_by_id,
            )
            if not ids:
                continue
            locus = None
            for fid in ids:
                locus = locus_by_id.get(fid)
                if locus:
                    break
            ck, cid, cs, ce = _build_composites(
                fr, s, e, types_by_id, tis_range=tis_range, flank_nt=flank_nt
            )
            out_rows.append(
                {
                    "orf_id": r.get("orf_id"),
                    "tran_id": str(tid),
                    "start_pos_tran": s,
                    "stop_pos_tran": e,
                    "locus_id": locus,
                    "feature_chain": ids,
                    "feature_types": types,
                    "slice_start": ss,
                    "slice_end": ee,
                    "event_pos": pos,
                    "comp_kind": ck,
                    "comp_feature_id": cid,
                    "comp_slice_start": cs,
                    "comp_slice_end": ce,
                }
            )

    pl.from_dicts(out_rows).write_parquet(out_parquet)
