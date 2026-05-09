# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ABO Streamlit Version** is a web application for blood group genetic analysis, specializing in ABO and RHD blood group genotyping from DNA sequencing data (AB1 chromatogram files and FASTA sequences). It performs sequence alignment, variant detection, and allele identification using BioPython and provides interactive visualization of results.

## Architecture

### Main Application
- **main.py**: Streamlit web application with three main analysis tabs:
  1. **Chromatogram Check for Heterozygotes**: Visual analysis of AB1 traces to detect heterozygous positions
  2. **Exon-based SNP Analysis**: Displays sequence alignments and detected variants from FASTA analysis
  3. **Allele Prediction**: Identifies possible alleles based on variant combinations using Cartesian product of exon results

### Core Analysis Modules (utils/)

#### Sequence Analysis
- **FASTA_analyzer.py**: Gap-aware pairwise alignment using BioPython's PairwiseAligner. Key features:
  - Multi-exon analysis of ABO sequences (exons 1-7)
  - Automatic exon detection and strand determination (forward/reverse)
  - Variant detection (SNPs, insertions, deletions)
  - Returns structured results with similarity scores for batch processing

- **ab1_analyzer.py**: AB1 file parsing (Sanger sequencing traces):
  - Extracts DNA traces (A, C, G, T channels) from AB1 binary files
  - Merges overlapping traces from multiple files
  - Detects heterozygous positions based on signal intensity ratio
  - Returns normalized trace data for chromatogram visualization

#### Allele & Phenotype Identification
- **abo_identifier.py**: Graph-based ABO allele identification:
  - Loads ISBT allele definitions as a network graph
  - Matches detected variants to known alleles
  - Supports heterozygous variant handling via IUPAC codes
  - Returns intersection of possible alleles across variants

- **rhd_analyzer.py**: RHD blood group genotyping:
  - Embeds WHO standard RHD1 (951bp) and RHD456 (~3336bp) reference sequences
  - Auto-detects amplicon region by alignment length and identity
  - Length-based decision for RHD1 (< 800bp = RhD+, indicating RHD gene present)
  - Identity-based decision for RHD456 (≥90% = RhD+)
  - Multi-amplicon support with consolidation logic for final verdict

- **isbt_handler.py**: ISBT database integration:
  - Loads reference ABO allele data from JSON
  - Suggests blood group phenotypes from detected variants
  - Maps variant combinations to ISBT phenotype recommendations

#### Blood Group System Registry
- **bloodgroup/registry.py**: Modular system for registering and accessing blood group systems
  - Extensible design for adding new blood group systems (RhCE, Kell, Duffy, Kidd, H)
  - Central point for system configuration and availability

#### Utilities
- **referece_loader.py**: Loads ABO reference data from JSON (gene_sequences.json)
  - Exon coordinates, CDS boundaries, reference sequences for all exons

### Data Files (utils/data/)
- **gene_sequences.json**: ABO reference sequences and exon metadata
- **abo_reference_ng006669.json**: NCBI NG_006669 reference with allele definitions
- **FUT1_enhanced_subgraph.gml, FUT2_enhanced_subgraph.gml, ABO_enhanced_subgraph.gml**: NetworkX graph files mapping variants to alleles
- **students.csv**: Team member information

## Development Setup

### Prerequisites
- Python 3.12+
- Virtual environment (recommended)

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
The app will open at `http://localhost:8501`

## Key Concepts & Patterns

### Exon-Based Analysis Strategy
ABO analysis works by:
1. Processing one exon per FASTA file (user uploads exon-specific FASTA files)
2. Aligning each to the corresponding ABO reference exon using local gap-aware alignment
3. Detecting variants (SNPs, insertions, deletions) in the aligned regions
4. Per-exon scoring: coverage (how much of reference is sequenced) and similarity (% identity)
5. Confirming exons with >90% similarity as "Confirmed"
6. Building Cartesian product of all confirmed exons to find consistent allele combinations

### Variant Handling
- IUPAC ambiguity codes (R=A/G, Y=C/T, etc.) represent heterozygous positions
- The system maps each detected variant to possible alleles via graph lookup
- Heterozygous variants expand the search space (e.g., "A_SNP" matches multiple allele nodes)
- Final alleles are the intersection of all variant-matched alleles

### RHD Logic
- RHD1 amplicon: short region (< 800bp), presence/absence determines RhD status
- RHD456 amplicon: longer region (> 2000bp), high identity (≥90%) = RhD+
- Multi-amplicon samples: consolidated verdict (all must agree on phenotype)
- Built-in WHO reference sequences ensure consistent results without external files

## Common Development Tasks

### Adding a New Blood Group System
1. Create analyzer module in utils/ (e.g., `kell_analyzer.py`)
2. Register system in `utils/bloodgroup/registry.py`
3. Add reference data to `utils/data/` as needed
4. Update main.py UI to include new system tabs/options

### Modifying FASTA Alignment Parameters
- Edit alignment parameters in `FASTA_analyzer.py.__init__()`:
  - `local_alignment_params` for exon alignment
  - `global_alignment_params` for full-sequence alignment
  - Scoring: match_score (positive), mismatch_score (negative), gap penalties

### Adjusting RHD Decision Thresholds
- In `rhd_analyzer.py`, update these constants:
  - `RHD1_MAX_LENGTH`: length threshold for RHD1/RHCE discrimination
  - `IDENTITY_RHD_POSITIVE`: identity % threshold for RhD+ call
  - `MAX_VARIANTS_RHD_POSITIVE`: variant count threshold for RhD+ status

### Debugging Sequence Alignment
- Set `DEBUG = True` in `FASTA_analyzer.py` to print alignment details
- The `_extract_aligned_sequences()` method parses BioPython alignment format
- Variant positions are identified by comparing aligned sequences character-by-character

## Dependencies
- **BioPython**: Sequence alignment, AB1 file parsing (PairwiseAligner, SeqIO)
- **Streamlit**: Web UI framework
- **Plotly**: Interactive chromatogram visualization
- **NumPy/Pandas**: Data manipulation
- **Matplotlib**: General plotting (backup for Plotly)
- **NetworkX**: Graph-based allele matching
- **SQLAlchemy**: (included, not actively used)

## Git Workflow
- Main branch (`main`) is the primary development branch
- Recent commits show focus on ABO+RHD integration and UI improvements
- Use conventional commit messages (e.g., "add RHD analysis", "fix FASTA logic")

## Known Issues & Notes
- AB1 merge logic handles overlapping regions by averaging signal intensity
- Heterozygote detection threshold (default 0.3 ratio) may need adjustment for weak secondary peaks
- FASTA files must be exon-specific (one exon per file) for accurate per-exon scoring
- RHD analysis expects ~2kb amplicons (RHD1 or RHD456 regions)
- Reference sequence files are embedded (no external database required)
