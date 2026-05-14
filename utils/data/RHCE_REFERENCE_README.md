# RHCE Reference Data

The RHCE analyzer requires a reference sequence to align patient reads against.
You must supply one of the files below before running RHCE analysis.

## Option A (preferred): GenBank file with exon map

Download **NG_009208** (RHCE genomic) from NCBI and save as:

```
utils/data/rhce_referance.gb
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NG_009208.1

This gives the analyzer an exon map so cDNA positions (c.48, c.178, c.676, etc.)
can be translated to genomic coordinates automatically.

## Option B (fallback): cDNA FASTA

Download **NM_020485** (RHCE mRNA) from NCBI and save as:

```
utils/data/rhce_cdna_reference.fasta
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NM_020485.7

With cDNA, the diagnostic SNP cDNA positions (c.48, c.178, c.676...) map
directly to string indices in the reference.

## Diagnostic SNPs covered (Asian-population focus)

| SNP        | Exon | Antigen call    | Notes                              |
|------------|------|-----------------|------------------------------------|
| c.48G>C    | 1    | C (vs c)        | Primary C/c marker (p.Ser16Cys)    |
| c.178A>C   | 2    | C (vs c)        | Confirmatory C marker              |
| c.203G>A   | 2    | C (vs c)        | Confirmatory C marker              |
| c.307C>T   | 2    | C (vs c)        | Confirmatory C marker (p.Leu103Ser)|
| c.676G>C   | 5    | E (vs e)        | Primary E/e marker (p.Ala226Pro)   |
| c.1025T>C  | 7    | Cat. EII partial| Asian-specific partial E variant   |
| c.1226A>G  | 8    | Cat. EIII partial| Asian-specific partial E variant  |

Reference allele = **c** and **e** (ISBT RHCE*01 = ce).

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology, Table 004 (RH).
- Daniels G. *Human Blood Groups* 3rd ed. — Chapter on RHCE molecular basis.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
