# E. coli Synteny Explorer

Comparative genomics workflow for exploring local protein synteny between *Escherichia coli* genomes.

## Overview

This project identifies a protein of interest in a reference genome, retrieves its immediate genomic neighbours, searches for homologous proteins in a target genome with BLASTP, and assembles the results into a synteny-oriented table.

It combines sequence parsing, genome annotation tables, command-line BLAST+ and Python-based data processing.

## Workflow

```text
Reference genome + protein of interest
                │
                ▼
       Identify neighbouring genes
                │
                ▼
       Extract protein sequences
                │
                ▼
          BLASTP against
          target genome DB
                │
                ▼
        Select best hits
                │
                ▼
 Compare genomic positions / strand
                │
                ▼
          Synteny table
```

## Technologies

- **Python**
- **Biopython** for FASTA and BLAST XML parsing
- **pandas** for genome annotation and result handling
- **BLASTP / BLAST+** for cross-genome protein homology searches
- **Jupyter + ipywidgets** for interactive exploration

## Analysis steps

For a selected protein, the workflow:

1. locates the protein from a `Protein_Id`, `Gene_Id` or `Gene_Name`;
2. determines its genomic position and strand;
3. retrieves the immediate upstream and downstream proteins;
4. extracts protein sequences;
5. runs BLASTP against a target-genome database;
6. parses BLAST XML output and identifies the best hit;
7. retrieves target-genome positions;
8. combines the information into a synteny-oriented table.

The implementation also accounts for gene strand when defining upstream and downstream neighbours.

## Repository structure

```text
.
├── projet.py
├── jupyter_notebook.ipynb
├── requirements.txt
├── .gitignore
└── README.md
```

## Data and BLAST databases

Genome annotation files, protein FASTA files and BLAST databases are not included in this repository. They should be generated from the relevant public *E. coli* genome resources before running the workflow.

The current implementation expects the following project layout:

```text
data/
├── Ecoli_genomes_refseq.xlsx
└── genomes/
    └── <assembly_accession>/
        ├── annotation_<assembly_accession>.tsv
        ├── protein.faa
        └── <assembly_accession>.*   # BLAST database files
```

BLAST+ must be installed locally and available to the workflow.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows:

```powershell
.venv\Scripts\activate
pip install -r requirements.txt
```

## Running the notebook

```bash
jupyter notebook
```

Open `jupyter_notebook.ipynb` and follow the interactive workflow from genome selection to synteny-table generation.

## Reproducibility

The repository documents the expected input structure and Python dependencies. The BLAST databases and genome files are intentionally kept outside the repository because they are generated from external genome resources and can be large.

The current code reflects the original implementation. A future refactor will separate the analysis into a reusable package, make the BLAST executable path configurable, add automated tests and provide a small redistributable example dataset.

## Skills demonstrated

**Bioinformatics:** comparative genomics · synteny · protein homology · genomic neighbourhoods  
**Programming:** Python · pandas · Biopython  
**Tools:** BLAST+ · Jupyter · ipywidgets  
**Good practice:** documented inputs · environment specification · explicit limitations

---

*Portfolio project — Bioinformatics & Computational Biology.*
