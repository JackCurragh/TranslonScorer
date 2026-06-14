# RDG-Flux Export Contract

The v1 RDG-Flux export is a position-level substrate table. It is not an ORF call table and does not collapse replicate identity.

Command:

```bash
python3 -m TranslonScorer.cli export-rdg-flux \
  --profiles <transcript_profiles.parquet> \
  --frame-support <frame_support.parquet> \
  --sample-id <sample_or_replicate_id> \
  --out <rdg_flux_frame_posteriors.parquet> \
  --annotation-source <annotation.gtf> \
  --transcriptome-fasta <transcriptome_or_genome.fa> \
  --psite-offset-model <offset_model_label_or_path>
```

If `--frame-support` is not supplied, the command can build frame support from the profile and one of `--cds`, `--annotation`, or `--annotation-dir`:

```bash
python3 -m TranslonScorer.cli export-rdg-flux \
  --profiles <transcript_profiles.parquet> \
  --annotation <annotation.gtf> \
  --sample-id <sample_or_replicate_id> \
  --out <rdg_flux_frame_posteriors.parquet> \
  --frame-method linear+hmm
```

Output columns:

- `sample_id`
- `transcript_id`
- `pos`
- `count`
- `p_frame0`
- `p_frame1`
- `p_frame2`
- `p_background`
- `frame_entropy`
- `effective_depth`
- `p_translated`
- `local_periodicity_score`

Coordinate convention:

- `pos` is the same 0-based transcript coordinate used by the transcript profile parquet.
- `transcript_id,pos` therefore joins directly to TranslonScorer profile outputs after renaming `tran_id` to `transcript_id`.

Probability convention:

- `p_frame0..2` are posterior mass for translated signal in each transcript frame.
- Frame-support `p0/p1/p2` are conditional frame probabilities. If the frame-support table includes `support_evidence`, export scales those conditional probabilities by that evidence.
- `p_translated = (1 - background_probability) * support_evidence` when `support_evidence` is present.
- If `support_evidence` is absent, `p_translated` falls back to `1 - background_probability` for codons with frame support and `0` for unsupported codons.
- `p_background` is the residual mass, `1 - p_translated`.
- `p_frame0 + p_frame1 + p_frame2 + p_background = 1` within floating point tolerance.
- Low-confidence frame evidence remains visible through high `frame_entropy`, low `local_periodicity_score`, and low `effective_depth`.

Metadata:

The sidecar JSON records coordinate system, P-site offset model, annotation source, FASTA, TranslonScorer version, model stage, normalization, sample ids, and frame method.

Partitioning:

Use `--partition-by-sample` to write Hive-style sample partitions:

```text
<out_dir>/
  metadata.json
  sample_id=<sample_a>/part.parquet
  sample_id=<sample_b>/part.parquet
```

Scaling:

- With precomputed Parquet frame support and single-file output, the exporter uses Polars lazy scan/sink and does not materialise the full RDG table as Python objects.
- If frame support is built on the fly from annotation, the frame-support build is still the heavier table-wise step. For genome-scale profiles, validate on a curated subset, then run full frame-support generation as a batch job and export from the precomputed frame-support parquet.
- Aggregated BigWig-derived transcript profiles work for v1 frame posterior export, but they only support global frame leakage correction because read length and read identity have already been collapsed.

Stage 3 extension:

Joint transcript-frame posterior export should be added as a separate sparse posterior layer with explicit residual mass. It should not replace this v1 per-position frame table.

First integration target:

```bash
python3 -m TranslonScorer.cli export-rdg-flux \
  --profiles /Users/jackt/projects/all-RiboSeq/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet \
  --annotation /Users/jackt/projects/all-RiboSeq/bench/hg38/gencode.v45.annotation.gtf \
  --sample-id ribocrypt_fwd_full \
  --out /Users/jackt/projects/all-RiboSeq/ribocrypt_fwd/ribocrypt_fwd_rdg_flux_frame_posteriors.parquet \
  --transcriptome-fasta /Users/jackt/projects/all-RiboSeq/bench/hg38/hg38.fa \
  --annotation-source /Users/jackt/projects/all-RiboSeq/bench/hg38/gencode.v45.annotation.gtf \
  --frame-method linear
```

The real profile currently has 157,586,237 rows, so the full export should be run as a batch job after frame-support generation has been validated on the curated locus subset.
