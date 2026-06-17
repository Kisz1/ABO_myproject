# Decision thresholds — provenance, rationale, and sensitivity

This document catalogs **every hard-coded numeric threshold** that gates a
call in the pipeline, why it has the value it does, and (where applicable) the
**measured** evidence that calls do not flip under reasonable variation.

Each threshold is justified in one of two ways, stated explicitly per row:

- **Cited** — fixed by an external standard or a definitional relationship
  (e.g. the Phred quality scale, ISBT allele nomenclature).
- **Empirical** — an operational gate tuned for this assay and validated by the
  sensitivity sweep in [`tools/threshold_sensitivity.py`](../tools/threshold_sensitivity.py)
  (results: [`threshold_sensitivity_report.md`](threshold_sensitivity_report.md)).
  These are *not* claimed to be uniquely correct; the evidence is that the
  call is **stable across a plateau** that comfortably contains the chosen value.

> **Live vs dormant.** The content-based router
> (`utils/bloodgroup/router.py`) is marked *SKETCH ONLY* and is **not wired
> into the live app** (the app routes by filename via `route_filename`). Its
> `MIN_ROUTE_IDENTITY` / `MIN_ROUTE_MARGIN` thresholds therefore gate no
> production call today; they are documented here for completeness and flagged
> as dormant. (The legacy dead RHD methods `determine_rhd_phenotype` and
> `_determine_vote`, which referenced undefined constants, have since been
> **removed**; the live RHD verdict path is `determine_rhd_phenotype_snp_based`
> + `_calculate_final_verdict` via `analyze_multiple_amplicons`.)

---

## A. Base quality (Phred) — **Cited**

| Constant | Value | File | Gates |
|---|---|---|---|
| `MIN_PHRED_AT_SNP` | 30 | `rhd_analyzer.py`, and identically in `rhce/kel/fy/jk/h/mns/di_analyzer.py` | Per-base quality required *at* a diagnostic SNP column before the SNP is trusted (AB1 inputs only). |
| `QTRIM_Q_THRESHOLD` | 20 | `file_processing.py` | Sliding-window end-trim / internal N-mask floor. |
| `QTRIM_WINDOW` | 10 | `file_processing.py` | Window length for the trim. |

**Why.** The Phred quality score is defined as *Q = −10·log₁₀(P_error)*
(Ewing & Green, *Genome Res.* 1998;8:186–194; Ewing, Hillier, Wendl & Green,
*Genome Res.* 1998;8:175–185). This makes the values definitional, not tuned:

- **Q30 ⇒ P_error = 10⁻³** (1 erroneous base call in 1000) — the standard
  "high-confidence base" floor used for variant calling. Applied at the single
  SNP column, where one wrong base changes the allele.
- **Q20 ⇒ P_error = 10⁻²** (1 in 100) — the conventional Sanger end-trim floor;
  used over a 10-bp window (not per base) so a single dip doesn't truncate a
  read, while a sustained low-quality run does.

---

## B. Heterozygote calling — **Empirical** (+ definitional structure)

| Constant | Value | File | Gates |
|---|---|---|---|
| `HET_RATIO_PIPELINE` | 0.30 | `file_processing.py` | Clinical secondary/primary peak ratio to call a 2nd allele. |
| `HET_RATIO_DEFAULT` | 0.10 | `ab1_analyzer.py` | Permissive analyzer default (overridden by the pipeline value). |
| `HET_HALF_WINDOW` | 3 | `ab1_analyzer.py` | ± trace samples scanned around each base index for the peak max. |
| `HET_MIN_SUM` | 0.8 | `ab1_analyzer.py` | Minimum summed 4-channel signal; skips baseline/gap regions. |
| `HET_MIN_MAJOR` | 0.5 | `ab1_analyzer.py` | Minimum major-peak height before a position is eligible. |

**Why.** A heterozygous Sanger position shows a secondary chromatogram peak; the
question is only *how tall* that peak must be relative to the major peak. This
is a sensitivity/specificity trade-off, not a quantity fixed by an external
standard, so it is justified **empirically**:

- The pipeline uses **0.30** (a secondary peak ≥30% of the major) — stricter
  than the analyzer default **0.10** — trading sensitivity for specificity so
  noise/bleed-through is not over-called as a heterozygous SNP in clinical
  output.
- The sensitivity sweep (§Sensitivity) drives a controlled synthetic trace
  across secondary-peak ratios 0.05–0.50 and confirms the gate behaves as a
  clean cutoff (called iff data ratio ≥ threshold), with the 0.30 value
  rejecting the low-ratio (≤0.25) noise band that 0.10 would accept.

`HET_MIN_SUM` / `HET_MIN_MAJOR` / `HET_HALF_WINDOW` are structural guards (a
real base must be present and clearly peaked before any ratio test), not
diagnostic cutoffs.

---

## C. RHD identity ladder — **Cited (allele definitions) + Empirical (gates)**

| Constant | Value | File | Gates |
|---|---|---|---|
| `IDENTITY_RHD_DELETION` | 60.0 | `rhd_analyzer.py` | Below ⇒ complete RHD gene deletion (European RhD−). |
| `IDENTITY_RHD_PARTIAL_D` | 85.0 | `rhd_analyzer.py` | Floor for partial-D SNP calls (DIVa, DVI) + borderline band. |
| `IDENTITY_RHD_WEAK_D` / `IDENTITY_RHD_POSITIVE` | 90.0 | `rhd_analyzer.py` | "Gene confidently present" floor for weak-D / DEL SNP calls. |
| `IDENTITY_RHD_STANDARD_D` | 95.0 | `rhd_analyzer.py` | Clean Standard D with zero diagnostic SNPs. |
| `MIN_AMPLICON_IDENTITY` | 75.0 | `rhd_analyzer.py` | Voting floor: amplicons below this vote *Inconclusive* (noise/off-target). |

**Why.** The *allele assignments* keyed off these tiers follow **ISBT 004
v6.4** RHD nomenclature (cited inline in `ISBT_RHD_ALLELES`), and the RHDψ
pseudogene logic follows **Singleton et al., *Blood* 2000;95:12–18** and **Chiu
et al. 2001**. The *identity gate values* themselves are operational (how much
intact RHD sequence to require before trusting a call) and are justified
**empirically** on the real n=9 panel: the sweep shows every patient verdict has
large headroom above its gate and the 9 calls stay 100% concordant across a wide
plateau around each production value (§Sensitivity). The `MIN_AMPLICON_IDENTITY`
= 75 floor sits above the ~60% incidental RHD/RHCE homology noted in code.

---

## D. Probe-based SNP calling (RHCE, KEL, FY, JK, H, MNS, DI) — **Empirical + structural**

These seven systems share an identical, already-named threshold set:

| Constant | Value | Gates | Justification |
|---|---|---|---|
| `MIN_LOCAL_IDENTITY_FOR_CALLING` | 85.0 | Local/probe window identity required to trust a SNP. | Empirical — same gene-present confidence level as the RHD partial-D floor; introns drag *whole-read* identity down, so calling is done in a local window. |
| `MIN_IDENTITY_FOR_CALLING` | 85.0 | Whole-read identity (reporting / `callable_read` flag). | Empirical — coarse; reporting only, the local gate is the real criterion. |
| `SNP_LOCAL_WINDOW_MIN_BASES` | 30 | Min aligned non-gap bases in the SNP window. | Structural — a local identity % is meaningless below a minimum denominator. |
| `SNP_PROBE_MIN_BASES` | 60 | Min aligned bases for a usable probe hit. | Structural — a short probe match can be coincidental. |
| `MIN_PHRED_AT_SNP` | 30 | Q-gate at the SNP column (AB1). | Cited — see §A. |

The 85% local-identity gate carries even more headroom than RHD here: a true
SNP-bearing read aligns to its probe at near-100% within the window regardless
of surrounding introns.

---

## E. ABO exon confirmation — **Empirical + structural**

| Constant | Value | File | Gates |
|---|---|---|---|
| `ABO_EXON_MIN_COVERAGE` | 80.0 | `FASTA_analyzer.py` | % of the exon reference spanned by the alignment. |
| `ABO_EXON_MIN_SIMILARITY` | 70.0 | `FASTA_analyzer.py` | % identity within the aligned region. |
| `ABO_ALN_MIN_SCORE` | 20 | `FASTA_analyzer.py` | Raw-score noise floor; below it coverage/similarity are zeroed. |

Two-part confirmation (span **and** identity) before variant calling: coverage
guards against a short read matching a fragment by luck; similarity guards
against a different exon aligning weakly. Both are operational gates.

---

## F. Content-based routing — **DORMANT (sketch only)**

| Constant | Value | File | Status |
|---|---|---|---|
| `MIN_ROUTE_IDENTITY` | 90.0 | `bloodgroup/router.py` | Not wired into live app; below ⇒ "unknown". |
| `MIN_ROUTE_MARGIN` | 5.0 | `bloodgroup/router.py` | Not wired into live app; top−runner-up below ⇒ "ambiguous". |

The margin gate exists because RHD and RHCE share ~90% identity; documented for
when/if the content router replaces filename routing.

---

## Sensitivity

See [`threshold_sensitivity_report.md`](threshold_sensitivity_report.md),
regenerated by:

```
python tools/threshold_sensitivity.py
```

It sweeps the RHD ladder and amplicon floor on the **real n=9 RHD panel** and
the heterozygote ratio on a **controlled synthetic chromatogram**, reporting
whether any call flips.

### Measured results (latest run)

**RHD panel — 9/9 concordant at production thresholds, with zero verdict flips
across every swept range:**

| Threshold (production value) | Range swept | Patients flipped anywhere in range |
|---|---|---|
| `IDENTITY_RHD_DELETION` (60.0) | 40 → 75 | none (9/9 throughout) |
| `IDENTITY_RHD_PARTIAL_D` (85.0) | 75 → 90 | none (9/9 throughout) |
| 90% gene-present tier (90.0) | 82 → 95 | none (9/9 throughout) |
| `IDENTITY_RHD_STANDARD_D` (95.0) | 90 → 99 | none (9/9 throughout) |
| `MIN_AMPLICON_IDENTITY` (75.0) | 60 → 85 | none (9/9 throughout) |

**Interpretation.** The final RhD+/− verdict is produced by multi-amplicon,
SNP-based *weighted voting* (`REGION_VOTE_WEIGHT`), not by a single identity
cutoff — so it is robust to the exact identity-gate value. Concretely, a
reviewer asking *"why 90%?"* can be answered with data: **the call does not
change anywhere between 82% and 95% on the real panel.** Per-patient identity
headroom is large (RhD+ patients reach 100%; the RhD− patient KT is called from
SNP/voting evidence, not a low identity).

**Heterozygote ratio (synthetic).** The gate is a clean minor/major cutoff
(call iff data ratio ≥ threshold). The production **0.30** rejects the ≤0.25
secondary-peak band that the permissive **0.10** default would accept —
the specificity gain the pipeline override buys.

**Honest limitation.** This panel is 8 clear RhD+ and 1 clear RhD−. It robustly
exercises the **RhD+/− verdict**, but it does **not** contain borderline
weak-D / partial-D / DEL serology cases that would stress the *inter-tier*
(85 / 90 / 95%) boundaries used for **allele sub-typing**. Those tier values are
independently anchored to ISBT 004 v6.4 SNP definitions (§C) and are exercised
by the synthetic-construct unit tests (`tests/test_rhd_*`). Extending the
sensitivity analysis to a weak-D/partial-D panel is the natural next validation
step.
