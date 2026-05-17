# FY (Duffy / ACKR1) Reference Data

The FY analyzer requires a reference sequence to align patient reads against.

**Important:** the GATA-1 promoter SNP (-67T>C, defines FY\*Null / Fy(a-b-))
sits 67 bp 5' of the ATG start codon and is therefore **not present in the
mRNA**. To call this clinically critical SNP you must use the genomic
reference (Option A). cDNA-only (Option B) will silently mark the GATA SNP
as uncoverable.

## Option A (REQUIRED for GATA SNP): Genomic GenBank

Download **NG_011626.3** (ACKR1 genomic, ~7 kb including the promoter) from
NCBI and save as:

```
utils/data/fy_referance.gb
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NG_011626.3

Or fetch via NCBI eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_011626.3&rettype=gb&retmode=text
```

With the genomic record the analyzer reads the `exon` and `CDS` features
to build an exon map and locate the ATG; the promoter SNP is then addressed
as `cds_start_genomic + (-67)`.

## Option B (fallback): cDNA FASTA

Download **NM_002036.4** (ACKR1 mRNA, ~1.4 kb) from NCBI and save as:

```
utils/data/fy_cdna_reference.fasta
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NM_002036.4

Or fetch via NCBI eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_002036.4&rettype=fasta&retmode=text
```

With cDNA the analyzer locates the ATG (`ATGGGCAACTGCCT` context) and
subtracts the 5' UTR offset. The c.125G>A and c.265C>T SNPs map directly.
**The GATA -67T>C SNP cannot be called** — it will appear in the per-SNP
table as `not covered` with reason "promoter SNP requires a genomic
reference".

## Diagnostic SNPs covered (Phase 1)

| SNP             | Position           | Antigen call      | Notes                                                                            |
|-----------------|--------------------|-------------------|----------------------------------------------------------------------------------|
| c.125G>A        | Exon 2             | FY\*A vs FY\*B    | Primary p.Gly42Asp discriminator                                                 |
| -67T>C (GATA)   | Promoter (-67/ATG) | FY\*Null silencer | Erythroid silencer — Fy(a-b-) when homozygous; genomic reference required        |
| c.265C>T        | Exon 2             | FY\*X (weak B)    | p.Arg89Cys — weak FYB expression modifier                                        |

Reference allele = FY\*A (`G` at c.125). Variant allele = FY\*B (`A` at c.125).

## Phenotype determination

Phenotype calls combine the three SNPs:

- **GATA hom alt** (T/T → C/C) → `Fy(a-b-)` regardless of A/B genotype
  (both haplotypes erythroid-silenced).
- **GATA hom ref** (or uncallable) → A/B genotype determines phenotype:
  - `AA` → `Fy(a+b-)`
  - `AB` → `Fy(a+b+)`
  - `BB` → `Fy(a-b+)`
- **GATA het + A/B het** → phase ambiguous; analyzer returns both options
  (`Fy(a+b-)` if GATA cis with FY\*B, `Fy(a-b+)` if GATA cis with FY\*A).
- **FY\*X present** (c.265 het or hom alt) → adds a "(weak)" tag to any
  Fy(b+) call.

## Reference-file audit

The FY analyzer assumes the reference deposit at c.125 is the **FY\*A base
`G`** (i.e. the cloned donor was Fy(a+b-) at this position). If a future
NCBI revision inherits a FY\*B clone (`A` at c.125), the regression test

```
tests/test_fy_analyzer.py::test_file_base_at_c125_matches_audit_note
```

will fail. Fix by swapping `ref_base`/`alt_base` and `ref_call`/`alt_call`
in `FY_DIAGNOSTIC_SNPS['c.125G>A']` and adding `parity_inverted: True`,
mirroring the RHCE c.178 treatment.

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology,
  Table 008 (FY / ACKR1).
- Tournamille C, Colin Y, Cartron JP, Le Van Kim C. *Disruption of a
  GATA motif in the Duffy gene promoter abolishes erythroid gene
  expression in Duffy-negative individuals.* Nat Genet. 1995;10(2):224-8.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
- UniProt Q16570 (ACKR1_HUMAN) — protein sequence for start-codon context.
