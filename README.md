# Automated Genomic Sequence Retrieval & Annotation Engine

An intelligent Python tool that parses natural query requests, fetches precise genomic sequences live via Ensembl REST APIs, and automatically annotates key sequence properties (GC content, CpG islands, and Open Reading Frames).

## Features
- **Natural Query Parsing:** Resolves human gene symbols and promoter window requests.
- **Automated REST Retrieval:** Directly queries Ensembl endpoints to generate standard FASTA files.
- **Analytical Annotation:** Calculates GC%, scans for CpG islands ($Obs/Exp > 0.6$), and predicts Open Reading Frames (ORFs).

## Usage
```bash
python3 Fetch.py
