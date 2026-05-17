# MNS (GYPA + GYPB) Reference Data

The MNS analyzer requires **two** reference sequences — GYPA carries the
M/N polymorphism, GYPB carries S/s — and both must be present.

## Required files

### GYPA — NG_007470.3

Download from NCBI and save as:

```
utils/data/mns_gypa_referance.gb
```

Direct link: https://www.ncbi.nlm.nih.gov/nuccore/NG_007470.3

eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_007470.3&rettype=gb&retmode=text
```

### GYPB — NG_007483.3

Save as:

```
utils/data/mns_gypb_referance.gb
```

Direct link: https://www.ncbi.nlm.nih.gov/nuccore/NG_007483.3

eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_007483.3&rettype=gb&retmode=text
```

## Diagnostic SNPs covered (Phase 1)

| SNP                | Gene  | Codon              | Effect                             | Notes                                                                  |
|--------------------|-------|--------------------|------------------------------------|------------------------------------------------------------------------|
| c.59C>T (GYPA)     | GYPA  | 20 (TCA→TTA)       | p.Ser20Leu                         | Primary **M / N** discriminator                                        |
| c.143C>T (GYPB)    | GYPB  | 48 (ACG→ATG)       | p.Thr48Met                         | Primary **S / s** discriminator — **parity-inverted** (deposit = s)    |

## Phenotype determination

Both axes combine into a 9-cell phenotype matrix:

| GYPA c.59 | GYPB c.143 (parity-inverted) | Phenotype |
|---|---|---|
| C/C (MM) | C/C (ss)  | M+N-S-s+ |
| C/C (MM) | C/T (Ss)  | M+N-S+s+ |
| C/C (MM) | T/T (SS)  | M+N-S+s- |
| C/T (MN) | C/C (ss)  | M+N+S-s+ |
| C/T (MN) | C/T (Ss)  | M+N+S+s+ — **phase-ambiguous** (MS/Ns or Ms/NS) |
| C/T (MN) | T/T (SS)  | M+N+S+s- |
| T/T (NN) | C/C (ss)  | M-N+S-s+ |
| T/T (NN) | C/T (Ss)  | M-N+S+s+ |
| T/T (NN) | T/T (SS)  | M-N+S+s- |

The **double-heterozygous (MN + Ss)** case has two haplotype options
(MS/Ns vs Ms/NS) that Sanger genotyping cannot phase — both are returned.

## Reference-file audits

**GYPA c.59:** NG_007470.3 carries the **M allele** (`C` at c.59, codon
TCA = Ser). No parity inversion needed.

**GYPB c.143:** NG_007483.3 carries the **s allele** (`C` at c.143,
codon ACG = Thr). The SNP entry is **parity-inverted**: `ref_call='s'`,
`alt_call='S'`, `parity_inverted: True`. If a future revision flips
to S (`T`), un-invert the entry; same treatment as JK c.838.

Regression tests in `tests/test_mns_analyzer.py` enforce both audits.

## Phase-2 roadmap (not in this module yet)

- **Codon-24 confirmatory marker** (c.71G>A + c.72T>A): the N allele
  changes codon 24 from `GGT` (Gly) to `GAA` (Glu) via *two
  simultaneous* base changes. c.71G>A alone gives `GAT` (Asp), which
  is neither M nor N — so it cannot be called as a single SNP. A
  haplotype-pair caller would address this.
- **U antigen detection** (mostly GYPB whole-gene deletion).
- **Miltenberger hybrid alleles** (Mi.III etc.) — requires GYPA-GYPB
  hybrid breakpoint detection, not point-mutation calling.

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology,
  Table 002 (MNS / GYPA / GYPB).
- Daniels G. *Human Blood Groups* 3rd ed. — MNS molecular basis chapter.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
- UniProt P02724 (GYPA), P06028 (GYPB) — protein sequences for codon mapping.
