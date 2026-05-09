# ABO + RHD Blood Group Analysis System

A web application for automated blood group genetic analysis from DNA sequencing data (AB1 chromatograms and FASTA sequences).

## Features

### ABO Analysis
- **Exon-based FASTA alignment** (exons 1-7)
- **Automatic exon detection** with strand determination
- **Variant detection** (SNPs, insertions, deletions)
- **Allele identification** using graph-based matching (ISBT alleles)
- **Phenotype prediction** (A, B, AB, O types)

### RHD Analysis
- **ISBT-informed allele classification** (RHD*01.01, RHD*01W.4, etc.)
- **Diagnostic SNP detection** (c.1227G>A, c.809T>G, c.1025T>C)
- **Multi-amplicon voting system** for reliable determination
- **RhD+ / RhD- determination** with allele subtyping
- **Weak D and DEL type identification** (especially common in East Asia)

## Quick Start

### Installation
```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### Running the Application
```bash
streamlit run main.py
```

Access at: http://localhost:8501

## Project Structure

```
ABO_streamlit_version/
├── main.py                 # Streamlit web application
├── config.py              # Configuration settings
├── requirements.txt       # Python dependencies
├── utils/                 # Core analysis modules
│   ├── FASTA_analyzer.py      # Sequence alignment & variant detection
│   ├── ab1_analyzer.py        # AB1 chromatogram parsing
│   ├── rhd_analyzer.py        # RHD genotyping (ISBT standard)
│   ├── abo_identifier.py      # ABO allele identification
│   ├── isbt_handler.py        # ISBT phenotype mapping
│   ├── reference_loader.py    # Reference sequence loading
│   └── data/              # Reference sequences & allele data
│       ├── gene_sequences.json
│       ├── abo_reference_ng006669.json
│       └── *.gml          # NetworkX allele graphs
├── sample_for_global_alignment/  # Test sample files
├── docs/                  # Documentation & reports
│   ├── READY_TO_USE.md
│   ├── REAL_PATIENT_TEST_REPORT.md
│   └── *.md              # Implementation notes
└── tests/                # Test scripts (development)
    ├── test_abo_analysis.py
    ├── test_rhd_phenotype.py
    └── *.py              # Other test files
```

## Usage

### 1. Upload Files

**For ABO Analysis:**
- Upload FASTA files (one per exon)
- Supported: exons 1-7

**For RHD Analysis:**
- Upload AB1 or FASTA files
- Single or multiple amplicons (RHD1, RHD456, etc.)

### 2. Analyze

Click the "Analyze" button in the sidebar

### 3. View Results

**Tab 1: Chromatogram** - AB1 trace visualization with heterozygote detection
**Tab 2: Exon Analysis** - Coverage, similarity, variants per exon
**Tab 3: Allele Prediction** - Possible ABO genotypes
**Tab 4: RHD Results** - RhD+ / RhD- determination with ISBT allele name

## Technical Details

### ABO Analysis
- Uses BioPython PairwiseAligner for local alignment
- Reference: NCBI NG_006669 ABO gene
- Outputs: Exon coverage, similarity %, detected variants, matching alleles

### RHD Analysis
- Uses ISBT allele nomenclature (ISBT 004 v6.4)
- Embedded WHO standard RHD1 (951bp) and RHD456 (3336bp) references
- SNP detection for major alleles: c.1227G>A, c.809T>G, c.1025T>C
- Multi-amplicon voting: combines results for final verdict
- Note: Complete ISBT classification available at https://blooddatabase.isbtweb.org/

### Decision Logic (RHD)
- **Identity ≥ 95% + no SNPs** → RhD+ (Standard D / RHD*01.01)
- **Identity ≥ 90% + c.1227G>A** → Weak D type 4 / DEL (RHD*01W.4)
- **Identity ≥ 90% + c.809T>G** → Weak D type 1 (RHD*01W.1)
- **Identity ≥ 90% + c.1025T>C** → Weak D type 2 (RHD*01W.2)
- **Identity < 85%** → RhD- (RHD deletion)

## Testing

Run test files in `/tests/` directory:
```bash
python tests/test_abo_analysis.py
python tests/test_rhd_phenotype.py
python tests/test_with_sample_files.py
```

See `/docs/` for detailed test reports.

## Limitations

### ABO Analysis
- Requires exon-specific FASTA files
- Designed for ~90%+ sequence quality
- Limited to standard ABO alleles (ISBT defined)

### RHD Analysis
- Simplified SNP detection (3 verified positions only)
- Does NOT detect: weak D edge cases, partial D variants, RHCE system, extended RHD variants
- Recommendations:
  - Use serological typing to confirm phenotype
  - For detailed variant analysis: contact reference laboratory
  - For complete allele classification: consult ISBT database

## Standards & References

- **ABO System:** ISBT ABO allele database
- **RHD System:** ISBT 004 RHD allele nomenclature (v6.4)
- **Alignment:** BioPython PairwiseAligner with WHO parameters
- **Sequencing:** Supports Sanger (AB1) and FASTA formats

## Contact & Support

For questions or issues:
1. Check `/docs/` for implementation details
2. Review test reports in `/docs/`
3. Consult ISBT standards at https://www.isbtweb.org/

---

**Version:** 1.0  
**Last Updated:** 2026-05-09  
**Status:** Production Ready
