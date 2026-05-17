# KEL Reference Data

The KEL analyzer requires a reference sequence to align patient reads against.
You must supply one of the files below before running KEL analysis.

## Option A (preferred): cDNA FASTA

Download **NM_000420.3** (KEL mRNA, ~2.5 kb) from NCBI and save as:

```
utils/data/kel_cdna_reference.fasta
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NM_000420.3

Or fetch via NCBI eutils (used by the bootstrap script):
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_000420.3&rettype=fasta&retmode=text
```

With cDNA, the diagnostic SNP cDNA positions map directly to string indices
in the reference (after the analyzer locates the ATG start codon and subtracts
the 5'UTR offset).

## Option B (optional): GenBank file with exon map

Download **NG_007492.1** (KEL genomic) from NCBI and save as:

```
utils/data/kel_referance.gb
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NG_007492.1

This gives the analyzer an exon map so cDNA positions can be translated to
genomic coordinates automatically. Not required when cDNA is present.

## Diagnostic SNPs covered (Phase 1)

| SNP        | Exon | Antigen call    | Notes                              |
|------------|------|-----------------|------------------------------------|
| c.578C>T   | 6    | K (vs k)        | Primary K/k marker (p.Thr193Met)   |

Reference allele = **k** (KEL\*02). Variant allele = **K** (KEL\*01).

### Future SNPs (Phase 2, drop into `KEL_DIAGNOSTIC_SNPS` when amplicons exist)

| SNP        | Exon | Antigen call    | Notes                              |
|------------|------|-----------------|------------------------------------|
| c.841C>T   | 8    | Kpa (vs Kpb)    | p.Arg281Trp                        |
| c.1790T>C  | 17   | Jsa (vs Jsb)    | p.Leu597Pro                        |

## Reference-file audit

`NM_000420.3` carries the **k-allele** base (`C`) at c.578 — i.e. the cloned
donor was K-negative, the common phenotype. No parity inversion is needed.

If NCBI ever revises this record and the deposited base flips to `T` (the
K-allele variant), the regression test
`tests/test_kel_analyzer.py::test_file_base_at_c578_matches_audit_note`
will fail. Fix by swapping `ref_base`/`alt_base` and `ref_call`/`alt_call`
in `KEL_DIAGNOSTIC_SNPS` and setting `parity_inverted: True`, mirroring the
RHCE c.178 treatment in `utils/rhce_analyzer.py`.

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology,
  Table 006 (KEL).
- Daniels G. *Human Blood Groups* 3rd ed. — Chapter on Kell molecular basis.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
- UniProt P23276 (KELL_HUMAN) — protein sequence for start-codon context.
