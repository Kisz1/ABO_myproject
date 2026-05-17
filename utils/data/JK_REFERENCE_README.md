# JK (Kidd / SLC14A1) Reference Data

The JK analyzer requires a reference sequence to align patient reads against.

**Important:** the Polynesian Jk(a-b-) splice-site SNP (`c.342-1G>A`, last
base of intron 5) is **not present in the mRNA** because the intron itself
is spliced out. To call this clinically critical SNP you must use the
genomic reference (Option A). cDNA-only (Option B) will silently mark
that SNP as uncoverable.

## Option A (REQUIRED for splice-site SNP): Genomic GenBank

Download **NG_011775.2** (SLC14A1 genomic, ~35 kb) from NCBI and save as:

```
utils/data/jk_referance.gb
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NG_011775.2

Or fetch via NCBI eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_011775.2&rettype=gb&retmode=text
```

With the genomic record the analyzer reads the `exon` and `CDS` features
to build an exon map and locate the ATG; intronic SNPs are then addressed
via `cDNA_anchor + intronic_offset`.

## Option B (fallback): cDNA FASTA

Download **NM_015865.7** (SLC14A1 mRNA, ~3 kb) from NCBI and save as:

```
utils/data/jk_cdna_reference.fasta
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NM_015865.7

With cDNA the analyzer locates the ATG (`ATGGAGGACAGCCCC` context for
SLC14A1 / UT-B) and subtracts the 5' UTR offset. The c.838 and c.871
SNPs map directly. **The c.342-1G>A intron-5 splice SNP cannot be called**
— it will appear in the per-SNP table as `not covered` with a clear reason.

## Diagnostic SNPs covered (Phase 1)

| SNP             | Position             | Antigen call         | Notes                                                                                |
|-----------------|----------------------|----------------------|--------------------------------------------------------------------------------------|
| c.838G>A        | Exon 9               | JK\*A vs JK\*B       | Primary p.Asp280Asn discriminator                                                    |
| c.342-1G>A      | Intron 5 (-1 bp)     | JK\*Null silencer    | Splice-acceptor abolished — Polynesian Jk(a-b-); genomic reference required          |
| c.871T>C        | Exon 9               | JK\*Null silencer    | p.Ser291Pro — Asian / Finnish Jk(a-b-) point mutation                                |

Reference allele = JK\*A (`G` at c.838). Variant allele = JK\*B (`A` at c.838).

## Phenotype determination

Phenotype calls combine the three SNPs:

- **Either null variant homozygous** → `Jk(a-b-)` regardless of A/B
  genotype (both haplotypes silenced).
- **No silencing detected** → A/B genotype determines phenotype:
  - `AA` → `Jk(a+b-)`
  - `AB` → `Jk(a+b+)`
  - `BB` → `Jk(a-b+)`
- **One null variant heterozygous + A/B het** → phase ambiguous; analyzer
  returns both options (`Jk(a+b-)` if null cis with JK\*B, `Jk(a-b+)` if
  null cis with JK\*A).
- **Compound heterozygous** (both nulls het in the same sample) is rare —
  flagged with LOW confidence in the combined-null call.

## Reference-file audit

The JK analyzer assumes the reference deposit at c.838 is the **JK\*A base
`G`** (the more common allele globally). If a future NCBI revision inherits
a JK\*B clone (`A` at c.838), the regression test

```
tests/test_jk_analyzer.py::test_file_base_at_c838_matches_audit_note
```

will fail. Fix by swapping `ref_base`/`alt_base` and `ref_call`/`alt_call`
in `JK_DIAGNOSTIC_SNPS['c.838G>A']` and adding `parity_inverted: True`,
mirroring the RHCE c.178 treatment.

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology,
  Table 009 (JK / SLC14A1).
- Irshaid NM, Henry SM, Olsson ML. *Genomic characterization of the kidd
  blood group gene: different molecular basis of the Jk(a-b-) phenotype
  in Polynesians and Finns.* Transfusion. 2000;40(1):69-74.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
- UniProt Q13336 (UT-B / SLC14A1) — protein sequence for start-codon context.
