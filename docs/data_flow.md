# Data flow

Two lanes that meet once. The annotation lane runs per `annotation_version` and
knows nothing about signal. The signal lane runs per cohort and knows nothing
about features. They meet at `_score_events_over_provider`, and only there.

```mermaid
flowchart TB
  subgraph A["① Annotation lane — once per annotation_version"]
    A1["sqlite · GTF · BED12 · bigBed · FASTA<br/><i>exactly one source</i>"]
    A5["context GTF<br/><i>optional</i>"]
    A2["feature_sources → blocks + translons"]
    A3["extract_events"]
    A4[("events/ · feature_event/ · event_overlap/")]
    A1 --> A2 --> A3 --> A4
    A5 -->|build_splice_context| A3
  end

  subgraph B["② Signal lane — three providers, one query contract"]
    subgraph B1["BAM"]
      BA["BAM(s)"] --> BB["BamSetProvider<br/><b>offsets calibrated at init</b><br/>per BAM × read-length, cached"]
    end
    subgraph B2["Matrix"]
      MA["partition dirs"] --> MB["Phase 1 · calibrate_offsets<br/>per sample × length + QC"]
      MB --> MC["Phase 2 · build-psite-index"]
      MC --> MD[("psite index — chrom-sharded<br/><b>offset baked into p_site</b>")]
      MD --> ME["MatrixProvider"]
      MA -.->|junctions + ledger<br/>read partition BAMs| ME
    end
    subgraph B3["Bigwig"]
      WA["stranded bigwig pairs<br/><b>offsets applied upstream,<br/>outside this tool</b>"] --> WB["BigwigSetProvider"]
    end
  end

  CON{{"coverage(regions, site, by_sample)<br/>size_factors()"}}
  BB --> CON
  ME --> CON
  WB --> CON

  A4 --> SC["_score_events_over_provider<br/>per chrom · merged padded windows"]
  CON --> SC
  MAP["mappability bigwig<br/><i>optional</i>"] -.->|diagnostic annotation only<br/>never gates a call| SC

  SC --> ASP["aspects · init | elongation | termination | junction"]
  ASP --> REC[("long-form score rows — _RECORD_SCHEMA")]
  REC --> CMP["_compose_per_translon"]
  CMP --> RPT[("per-translon report")]
```

## The three providers are not interchangeable

Every provider satisfies `coverage()` + `size_factors()`. Beyond that they differ,
and the differences are load-bearing rather than incidental.

| | BAM | Matrix | Bigwig |
|---|---|---|---|
| `coverage()` | ✅ | ✅ | ✅ |
| `junction_support()` | ✅ | ✅ | ❌ no CIGAR |
| `mappability_ledger()` | ✅ | ✅ | ❌ |
| **where offsets are applied** | provider init | index build | **upstream, outside the tool** |
| per-sample resolution | ✅ | ✅ | per-file only |

**Bigwig cannot score junction events at all.** `_warn_bigwig_cannot_score_junctions`
exists for this. Those events are unscorable, not unsupported, and the distinction
survives into the output.

**Offsets enter at three different points, and only two are inside the tool.**
The BAM path calibrates per BAM × read-length at provider init. The matrix path
bakes the offset into `p_site` at index build, which is why `MatrixProvider`
hard-errors without `psite_index_dir` — a flat offset across read lengths pins
`elong_in_frame` to the ~0.33 random floor. The bigwig path inherits whatever was
done upstream and cannot verify it.

Frame-sensitive metrics are therefore **not comparable across providers**, and a
threshold calibrated on one path is invalid on another.

## Invariants

1. An event is scored **once** per (provider, thresholds_version). Every feature
   claiming that event reads the same row.
2. `extract_events` is the only step that reads a feature source. Scoring never
   re-reads annotation.
3. The mappability track annotates; it never gates.
4. `eligibility` ("could we test this") is never conflated with `call`
   ("did it pass").
5. Thresholds are a versioned input. Changing one re-derives calls without
   re-measuring coverage.
