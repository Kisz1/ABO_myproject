# Methodology

## 1. Overview

We developed a web-based platform for **molecular blood-group genotyping from
Sanger sequencing data**. The system accepts raw chromatogram traces (`.ab1`)
and assembled sequences (`.fasta`) and reports predicted phenotypes and ISBT
alleles for nine blood-group systems: **ABO, RHD, RHCE, Kell (KEL), Duffy (FY),
Kidd (JK), H (FUT1), MNS (GYPA/GYPB), and Diego (SLC4A1)**.

The pipeline is implemented in Python with a Streamlit front end. Sequence
alignment throughout uses BioPython's `PairwiseAligner`. All variant
coordinates and phenotype calls are anchored to NCBI RefSeq references and
interpreted against the **ISBT (International Society of Blood Transfusion)**
allele tables.

The workflow has four stages: (i) input acquisition and signal processing,
(ii) reference alignment and variant calling, (iii) per-system phenotype/allele
interpretation, and (iv) consolidation and reporting.

---

## 2. Input Acquisition and Signal Processing

### 2.1 Two input modalities

Each sample can be supplied as either a **FASTA** read (already base-called) or
a raw **AB1 chromatogram**. FASTA reads enter the alignment stage directly. AB1
files are first processed to recover the base-called sequence, the four-channel
dye-signal traces, and (where required) per-base quality scores.

### 2.2 AB1 chromatogram parsing

For each AB1 file the four dye channels (A/C/G/T) are extracted from the raw
`DATA9–DATA12` arrays and mapped to nucleotides using the file's `FWO_` channel
order, with peak locations taken from `PLOC2`. Two read paths are used:

- **Orientation-aware path (ABO):** the called sequence is aligned to the ABO
  reference in both orientations, and the better-scoring strand is selected;
  reverse reads are reverse-complemented, with signal arrays and peak indices
  remapped consistently.
- **Quality-preserving path (RHD and extended systems):** the raw trace is read
  together with its Phred quality vector, and strand selection is deferred to
  the downstream analyzer (which aligns bidirectionally).

### 2.3 Signal normalization and heterozygote detection

Each channel is normalized per-channel by its 99th-percentile intensity to make
peak heights comparable. **Heterozygous positions** are detected directly from
the chromatogram: at each called peak we compare the major and minor channel
intensities within a small window (±3 samples), and flag a site as heterozygous
when the minor/major ratio exceeds an adjustable threshold (default 0.30),
subject to minimum total-signal and minimum major-peak gates to suppress noise.

### 2.4 Quality trimming and IUPAC encoding (RHD path)

For the quality-preserving path, reads undergo **sliding-window Phred trimming**:
5′ and 3′ ends are trimmed until a window of bases all meet the quality
threshold (default Q20, window 10), and internal low-quality bases are masked to
`N`. Detected heterozygous peaks are then re-encoded as **IUPAC ambiguity codes**
(e.g. A/G → R) in the trimmed sequence, so that heterozygous diagnostic SNPs are
preserved for downstream calling rather than being collapsed to a single base.

---

## 3. Reference Alignment and Variant Calling

### 3.1 Gap-aware pairwise alignment

All comparisons use BioPython's `PairwiseAligner` with explicit, consistent
scoring (match +2, mismatch −0.5 to −1, gap open −2 to −5, gap extend −0.5).
Local (Smith–Waterman) alignment is used for short amplicons and probe regions;
global alignment is used for long references. Using a **gap-aware aligner rather
than position-by-position comparison** eliminates the coordinate frame-shifts
and false-SNP cascades that arise when insertions or deletions are present.

### 3.2 ABO: exon-based variant calling

ABO reads are aligned exon-by-exon against the reference exon sequences. For
each exon we compute coverage and identity, and an exon is **confirmed** only
when it passes a two-step gate (coverage ≥ 80%, identity ≥ 70%). Variants are
then extracted from the confirmed alignment and classified as **SNP, insertion,
or deletion**, each carrying genomic and ISBT (c.) coordinates derived from the
exon's CDS offset. Both strands are evaluated and the winning strand is selected
per exon by alignment similarity.

### 3.3 ABO: allele identification

Confirmed variants (with IUPAC ambiguity expanded into candidate bases) are
matched against an allele/variant knowledge base. For each variant the set of
compatible alleles is retrieved, and the **intersection across all observed
variants** yields the candidate allele set. Variants that match no known allele
node are reported separately as *unknown* variants for expert review.

---

## 4. RHD Analysis

RHD is treated separately because RhD status depends on gene **presence/absence**
and on a panel of clinically important weak-D / partial-D variants, and because
RHD shares ~90% identity with RHCE (risking cross-amplification).

### 4.1 Region-aware identity classification

References for four amplicon regions (RHD1, RHD456, RHD7, RHD9) are derived as
verified subsequences of the RefSeqGene **NG_007494.1 / LRG_796**. Each read is
aligned to all regions and assigned to the **best-fitting region by alignment
identity** rather than by length, which is robust to the RHCE homology.

### 4.2 Diagnostic-SNP-based phenotype calling

Within the assigned region the analyzer combines overall identity with detection
of region-specific **diagnostic SNPs**, each annotated with zygosity:

- identity < 60% → complete RHD deletion (European RhD−);
- identity ≥ 90% with no diagnostic SNP → standard D antigen (RhD+);
- specific SNPs map to named clinical variants — e.g. c.1227G>A (Weak D type 4 /
  Asian DEL), c.809T>G (Weak D type 1), c.1025T>C (Weak D type 2), c.1154G>C
  (Weak D type 3), c.602C>G (Partial D IVa), c.667T>G (Partial D VI);
- borderline identities (85–90%) are flagged as variant/inconclusive.

An RHD-Ψ (pseudogene) screen is run before identity-based calling.

### 4.3 Multi-amplicon voting

When several RHD amplicons are provided, per-amplicon calls are consolidated
into a single verdict: unanimous positive or negative calls are reported as
*confirmed*, while disagreement is reported as *inconclusive — confirmation
required*. Region weights and an identity floor prevent low-quality amplicons
from polluting the vote.

---

## 5. Extended Systems (RHCE, KEL, FY, JK, H, MNS, DI)

The seven extended systems share a common **probe-based SNP-calling** design.

### 5.1 Probe-based, intron-tolerant SNP calling

For each ISBT diagnostic SNP, a short reference probe (≈ ±50 bp) centred on the
SNP is locally aligned against the read. Because each SNP is called from its own
local probe, the method is **intron-tolerant**: genomic Sanger amplicons can be
called against cDNA references, and promoter SNPs (e.g. the Duffy GATA −67T>C
silencing variant) can be called against a genomic reference using the same code
path. The base at the SNP column is read off the alignment, decoded through
IUPAC, and gated by local identity and, where AB1 quality is available, a
per-base Phred floor (Q30).

### 5.2 Zygosity and antigen-axis calling

Each SNP is called as homozygous-reference, heterozygous, or homozygous-alternate,
then mapped to antigen axes defined by ISBT (e.g. FY c.125G>A → Fy^a/Fy^b;
KEL c.578C>T → K/k; MNS GYPA c.59 → M/N and GYPB c.143 → S/s). Reference-allele
"parity" is explicitly tracked so that, if an NCBI reference deposit carries the
alternate allele, calls remain correct.

### 5.3 Multi-read consensus voting

When multiple reads cover the same system, each SNP is called per read and a
**consensus vote** is taken across reads. Agreement level sets a confidence label
(HIGH/MEDIUM/LOW), the number of callable vs. total reads is tracked, and SNPs
with no usable read are reported as *no-call*. Where two SNPs cannot be phased
from independent reads, the result is flagged **phase-ambiguous** and the
possible ISBT haplotype pairs are enumerated rather than forcing a single call.

---

## 6. Unified Auto-Routing

To remove manual per-system file assignment, an optional **auto-detect uploader**
routes each read to the correct system by two complementary mechanisms:

1. **Filename routing (fast path):** filenames are matched against lab
   amplicon-naming conventions (RHCE, RHD, ABO, KEL56, FY12, JK78, MIA234,
   DI1819, etc.), with specific patterns ordered before generic ones.
2. **Sequence-content routing:** each read is locally aligned (both strands)
   against every system's reference and ranked by identity. A read is *routed*
   when the top identity ≥ 90% **and** its margin over the runner-up ≥ 5 pp;
   a small margin (e.g. the RHD/RHCE homology case) is flagged *ambiguous* for
   manual override, and a low top identity is flagged *unknown*.

Routed reads are appended to the existing per-system input lists, so the
downstream analyzers are unchanged. The manual per-system uploaders are retained
as an override path.

---

## 7. Phenotype Interpretation and Reporting

Per-system results are passed through an **ISBT data handler** that suggests
phenotypes and matching alleles from the observed variant set. Results from all
nine systems are then consolidated into a single normalized summary, each system
carrying its status (analyzed / missing-data / inconclusive), predicted
phenotype, candidate alleles, confidence, and read-support counts. The interface
renders a consolidated final blood-group result plus per-system detail panels
(per-SNP consensus, alignment views, and per-read breakdowns).

---

## 8. Implementation

The platform is implemented in **Python 3.12** using **Streamlit** (interface),
**BioPython** (AB1 parsing and pairwise alignment), **NumPy/SciPy** (signal
processing and peak handling), and **Plotly/Matplotlib** (interactive
chromatogram and alignment visualization). Each analyzer is independently unit-
tested (including tests against real RHD patient samples and reference audits
that guard the reference-allele parity assumptions).
