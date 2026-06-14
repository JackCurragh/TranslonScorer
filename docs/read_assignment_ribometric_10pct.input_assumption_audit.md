# Read-Assignment Input Assumption Audit: read_assignment_ribometric_10pct

![Input assumption audit](../notebooks/read_assignment_ribometric_10pct.input_assumption_audit.png)

## Summary

- Alignments with at least one tested offset inside annotated CDS: 564,996.
- Best global offset from alignment start: 13.
- Best global annotated-CDS frame fraction: 0.368.
- Offset-0 annotated-CDS frame fraction: 0.365.
- Gate threshold: 0.550.
- Decision: this substrate **fails** the frame-validation input gate.

## Interpretation

This audit asks whether the coordinate used by the read-assignment prototype behaves like a P-site coordinate. It scans offsets from the BAM alignment start and measures how concentrated annotated-CDS positions are in one reading frame.

The substrate does not have strong enough global frame concentration for biological validation of frame-aware assignment. It can still test candidate-table construction, EM plumbing, and ambiguity reporting, but posterior movements from this substrate should not be treated as evidence that frame resolved transcript origin.

## Top Global Offsets

| offset | frame0 | frame1 | frame2 | max_frame | best_frame | CDS positions |
|---:|---:|---:|---:|---:|---:|---:|
| 13 | 0.368 | 0.342 | 0.290 | 0.368 | 0 | 561487 |
| 16 | 0.368 | 0.342 | 0.290 | 0.368 | 0 | 562652 |
| 19 | 0.367 | 0.342 | 0.291 | 0.367 | 0 | 563815 |
| 14 | 0.291 | 0.367 | 0.342 | 0.367 | 1 | 561971 |
| 10 | 0.367 | 0.342 | 0.291 | 0.367 | 0 | 557939 |
| 17 | 0.291 | 0.367 | 0.342 | 0.367 | 1 | 563069 |
| 20 | 0.291 | 0.367 | 0.342 | 0.367 | 1 | 564219 |
| 15 | 0.342 | 0.291 | 0.367 | 0.367 | 2 | 562321 |
| 7 | 0.367 | 0.342 | 0.291 | 0.367 | 0 | 555605 |
| 18 | 0.342 | 0.291 | 0.367 | 0.367 | 2 | 563450 |
| 11 | 0.292 | 0.367 | 0.341 | 0.367 | 1 | 558609 |
| 21 | 0.342 | 0.291 | 0.367 | 0.367 | 2 | 564648 |

## Best Offsets For Common Read Lengths

| read_length | best_offset | max_frame | best_frame | alignments |
|---:|---:|---:|---:|---:|
| 31 | 13 | 0.369 | 0 | 120965 |
| 30 | 12 | 0.354 | 0 | 120762 |
| 29 | 10 | 0.364 | 0 | 95076 |
| 32 | 13 | 0.390 | 0 | 89885 |
| 33 | 13 | 0.398 | 0 | 48750 |
| 28 | 13 | 0.356 | 0 | 43871 |
| 34 | 13 | 0.412 | 0 | 15890 |
| 27 | 9 | 0.368 | 0 | 13365 |
| 26 | 6 | 0.355 | 0 | 5232 |
| 35 | 10 | 0.409 | 0 | 2382 |
| 25 | 0 | 0.348 | 0 | 2143 |
| 24 | 0 | 0.390 | 0 | 1154 |
| 17 | 0 | 0.354 | 0 | 1047 |
| 18 | 2 | 0.366 | 0 | 900 |
| 16 | 15 | 0.350 | 0 | 540 |
