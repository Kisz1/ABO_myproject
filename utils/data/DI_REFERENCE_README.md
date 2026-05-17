# Diego (DI / SLC4A1) Reference Data

The Diego analyzer requires the SLC4A1 reference sequence.

## Option A (Recommended): Genomic GenBank

Download **NG_007498.1** (SLC4A1 RefSeqGene LRG_803, ~27 kb) from NCBI
and save as:

```
utils/data/di_referance.gb
```

Direct link: https://www.ncbi.nlm.nih.gov/nuccore/NG_007498.1

eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_007498.1&rettype=gb&retmode=text
```

## Option B (fallback): cDNA FASTA

Download **NM_000342.4** (SLC4A1 mRNA, ~3 kb) from NCBI and save as:

```
utils/data/di_cdna_reference.fasta
```

Direct link: https://www.ncbi.nlm.nih.gov/nuccore/NM_000342.4

The c.2561 SNP is in exon 19 of the coding sequence, so cDNA works
fine — but the genomic record is preferred for consistency and to
support future Phase-2 SNPs in other exons.

## Diagnostic SNPs covered (Phase 1)

| SNP        | Codon              | Effect                | Notes                                                               |
|------------|--------------------|-----------------------|---------------------------------------------------------------------|
| c.2561T>C  | 854 (CTG↔CCG)      | p.Pro854Leu           | Primary **Di(a) / Di(b)** discriminator — **parity-inverted entry** |

The lab amplicon **DI1819** (SLC4A1 exons 18-19) covers this position
and is the standard genotyping target.

## Phenotype determination

Single-axis system, three genotypes:

| c.2561 (parity-inverted) | Genotype | Phenotype |
|---|---|---|
| C/C | bb | Di(a-b+) |
| C/T | ab | Di(a+b+) |
| T/T | aa | Di(a+b-) |

## Reference-file audit

NG_007498.1 deposits the **Di-b allele** (`C` at c.2561, codon 854 =
CCG = Pro). This is the more common allele globally (>90%). The SNP
entry in `DI_DIAGNOSTIC_SNPS` is **parity-inverted**: `ref_call='b'`,
`alt_call='a'`, `parity_inverted: True`. If a future revision flips
to a Di-a clone (`T` at c.2561), un-invert the entry — same fix
pattern as JK c.838 and MNS GYPB c.143.

The regression test
`tests/test_di_analyzer.py::test_file_base_at_c2561_matches_audit_note`
enforces this audit.

## Phase-2 roadmap (not in this module yet)

- **Wr(a)/Wr(b)** at `c.1972G>A` (p.Glu658Lys) — exon 16 of SLC4A1.
  Not covered by DI1819 amplicons; would need a separate Wright-region
  amplicon and a second SNP entry. Wr(a) is extremely rare; Wr(b) is
  near-universal.
- Other Diego antigens (ELO, Wd(a), Rb(a), WARR) — all extremely rare
  population-specific variants; out of scope for routine genotyping.

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology,
  Table 010 (DI / SLC4A1).
- Daniels G. *Human Blood Groups* 3rd ed. — Diego molecular basis chapter.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
- UniProt P02730 (B3AT_HUMAN / SLC4A1) — protein sequence for codon mapping.
