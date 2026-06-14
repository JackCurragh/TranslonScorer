# Profile Disambiguation Walkthrough: read_assignment_ribometric_10pct

This report accompanies `notebooks/read_assignment_profile_disambiguation_walkthrough.ipynb`.

The plots compare profiles from selected ambiguous reads before and after disambiguation.

Important caveat: the RiboMetric BAM fails the current P-site/frame-signal input gate, so these are software and ambiguity diagnostics rather than biological validation plots.

Summary table: `../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.case_method_totals.csv`

## Plots

### LINC01783 / ENSG00000228549

Large movement in a tight lncRNA annotation block. Useful for ambiguity reporting, not coding-frame validation.

![LINC01783 / ENSG00000228549](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.LINC01783___ENSG00000228549.png)

### CIC

Coding-associated, but the reads span many loci. This mainly shows mappability ambiguity.

![CIC](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.CIC.png)

### ESRRAP2

Pseudogene or parent-gene ambiguity; the profile movement should not be counted as isoform resolution.

![ESRRAP2](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.ESRRAP2.png)

### OLFM1

Coding-associated multi-locus ambiguity with several plausible gene labels.

![OLFM1](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.OLFM1.png)

### PTCH1

Same-gene ambiguity, but the reviewed positions lack local independent frame evidence.

![PTCH1](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.PTCH1.png)

### RCE1

Same-gene ambiguity with no matched local frame support in this probe.

![RCE1](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.RCE1.png)

### H1-2 / H1-4

A compact paralog ambiguity; useful for locus-origin reporting, not transcript isoform validation.

![H1-2 / H1-4](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.H1-2___H1-4.png)

### NDUFS5

The right shape for a same-gene isoform case, but underpowered in this probe.

![NDUFS5](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.NDUFS5.png)

### POGK

Strong local frame support, but the read is cross-gene compatible.

![POGK](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.POGK.png)

### OAZ2 / IL17D

Single-read cross-gene ambiguity; visually useful, not a validation win.

![OAZ2 / IL17D](../notebooks/read_assignment_ribometric_10pct.profile_disambiguation.OAZ2___IL17D.png)
