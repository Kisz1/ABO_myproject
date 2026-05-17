# H (FUT1 / Bombay) Reference Data

The H analyzer requires a reference sequence to align patient reads against.

**Clinical context:** Loss-of-function variants on both FUT1 haplotypes
produce the **Bombay phenotype (Oh)** — RBCs lack H, A, and B antigens.
Bombay recipients **reject all conventionally ABO-typed donor blood** and
require Bombay-compatible units. Detecting Bombay carriers (Para-Bombay)
is therefore clinically critical pre-transfusion.

## Option A (Recommended): Genomic GenBank

Download **NG_007510.2** (FUT1 RefSeqGene, ~14 kb; bundles FGF21 + FUT1 +
IZUMO1) from NCBI and save as:

```
utils/data/h_referance.gb
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NG_007510.2

Or fetch via NCBI eutils:
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_007510.2&rettype=gb&retmode=text
```

The analyzer filters features by gene name (`FUT1` / `H` / `FT1` / `BCA`)
when loading, so the neighbouring FGF21 and IZUMO1 features are correctly
ignored. Without that filter, the first CDS encountered (FGF21 on the
minus strand) would become `cds_start_genomic` and every cDNA SNP position
would be mis-mapped.

## Option B (fallback): cDNA FASTA

Download **NM_000148.4** (FUT1 mRNA, ~3 kb) from NCBI and save as:

```
utils/data/h_cdna_reference.fasta
```

Direct link:
https://www.ncbi.nlm.nih.gov/nuccore/NM_000148.4

With cDNA the analyzer locates the ATG (`ATGTGGCTCCGGAGC` context for
FUT1; protein begins M-W-L-R-S per UniProt P19526) and subtracts the
5' UTR offset.

## Diagnostic SNPs covered (Phase 1)

All three are point mutations in the single FUT1 coding exon. Each was
verified against NG_007510.2 codons before being added to the analyzer.

| SNP        | Codon        | Effect                                     | Severity | Significance                                                                          |
|------------|--------------|--------------------------------------------|----------|---------------------------------------------------------------------------------------|
| c.725T>G   | 242 (CTC→CGC) | p.Leu242Arg                                | strong   | **Indian Bombay** — most common FUT1 loss-of-function variant worldwide               |
| c.586C>T   | 196 (CAG→TAG) | p.Gln196\* (nonsense)                      | strong   | Truncates FUT1 → classic Bombay null                                                  |
| c.460T>C   | 154 (TAC→CAC) | p.Tyr154His                                | weak     | h2 weak allele (FUT1\*01M.01) — reduced activity, Para-Bombay when homozygous         |

## Phenotype determination

The H system is binary at its core (H+ vs H-), modulated by partial
("weak") alleles. The analyzer aggregates all three SNPs into one of
three states:

- **H+** — all callable null SNPs are hom_ref → normal H expression.
- **Bombay (Oh) — H-negative** — any **strong** null SNP is homozygous →
  both haplotypes silenced; **rejects standard ABO blood**.
- **Para-Bombay** — either:
  - the weak allele (c.460) is homozygous → reduced but residual activity, **OR**
  - any null SNP is heterozygous → one haplotype affected, other functional.

## Reference-file audit

NG_007510.2 carries the wild-type FUT1 reference base at all three SNP
sites (T at c.725, C at c.586, T at c.460). No parity inversion is
needed. If a future NCBI revision flips any of these, regression tests
in `tests/test_h_analyzer.py` will fail with a clear message.

## Notes on Phase 1 scope

This module only handles **point mutations**. Several clinically
important Bombay alleles are **deletions** (c.547delC, c.880-882delTTC,
etc.) and are not yet supported — they would require an indel-aware
analyzer extension. Para-Bombay is also influenced by FUT2 / Lewis-type
secretor status (a separate gene), which is out of scope for the H
analyzer.

## Reference sources

- ISBT Working Party on RBC Immunogenetics & Blood Group Terminology,
  Table 018 (H / FUT1).
- Storry JR, Olsson ML. *The ABO blood group system revisited: a review
  and update.* Immunohematology. 2009;25(2):48-59.
- Kelly RJ, Ernst LK, Larsen RD, et al. *Molecular basis for H blood
  group deficiency in Bombay (Oh) and para-Bombay individuals.* Proc Natl
  Acad Sci U S A. 1994;91(13):5843-7.
- ISBT Allele Database: https://blooddatabase.isbtweb.org/
- UniProt P19526 (FUT1\_HUMAN) — protein sequence for start-codon context.
