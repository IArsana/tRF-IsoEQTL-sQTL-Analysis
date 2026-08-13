# PCa-tRFQTL

### Integrative Analysis of Prostate Cancer Risk Variants, tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

[![Status](https://img.shields.io/badge/status-in%20development-orange)]()
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

A reproducible bioinformatics pipeline for integrating prostate cancer genome-wide association study (GWAS), tRNA-derived fragment quantitative trait loci (tRF-QTL), isoform expression QTL (iso-eQTL), and splicing QTL (sQTL) data to identify shared regulatory loci associated with prostate cancer risk.

---

## Overview

Prostate cancer (PCa) is a genetically complex disease in which inherited genetic variation can influence disease susceptibility through diverse regulatory mechanisms.

While genome-wide association studies (GWAS) have identified numerous prostate cancer susceptibility loci, the functional mechanisms linking these variants to disease phenotypes remain incompletely understood.

Recent studies have demonstrated that genetic variants can regulate:

- tRNA-derived fragment (tRF) expression,
- transcript isoform expression,
- alternative splicing,
- and other layers of post-transcriptional regulation.

However, the potential relationship between **tRF genetic regulation and transcript isoform/splicing regulation at prostate cancer risk loci** remains relatively unexplored.

This project aims to investigate this regulatory landscape through an integrative QTL framework.

---

## Research Question

> **Do prostate cancer risk variants regulate tRNA-derived fragments and transcript isoform/splicing patterns through shared genetic regulatory mechanisms?**

The central hypothesis is that a subset of prostate cancer susceptibility variants may act as regulatory variants affecting both tRF expression and alternative transcript regulation.

Conceptually:

```text
                                         Prostate Cancer GWAS
                                                  │
                                                  ▼
                                         PCa risk variants
                                                  │
                                    ┌─────────────┴─────────────┐
                                    │                           │
                                    ▼                           ▼
                                 tRF-QTL                  iso-eQTL / sQTL
                                    │                           │
                                    ▼                           ▼
                                   tRF                  Transcript isoform /
                                                        alternative splicing
                                    │                           │
                                    └─────────────┬─────────────┘
                                                  ▼
                                        Shared regulatory loci
                                                  │
                                                  ▼
                                        Functional interpretation
                                                  │
                                                  ▼
                                         Clinical validation
```

## Objectives
### Primary Objective

To identify prostate cancer-associated genetic loci that simultaneously influence tRNA-derived fragment expression and transcript isoform/splicing regulation.

### Secondary Objectives
1. Identify prostate cancer risk variants associated with tRF expression.
2. Identify prostate cancer risk variants associated with transcript isoform expression.
3. Identify prostate cancer risk variants associated with alternative splicing.
4. Determine overlap between tRF-QTL and iso-eQTL/sQTL signals.
5. Investigate linkage disequilibrium (LD) relationships between associated variants.
6. Evaluate whether overlapping signals are consistent with shared causal variants using colocalization analysis.
7. Identify candidate genes and biological pathways associated with prioritized loci.
8. Evaluate the clinical relevance of prioritized tRFs and genes using prostate cancer datasets.
9. Validate candidate regulatory associations using independent datasets where available.

## Study Design

The analysis is organized into several major stages:

```text
                                    ┌─────────────────────────────┐
                                    │       Raw Data Sources      │
                                    ├─────────────────────────────┤
                                    │ Prostate Cancer GWAS        │
                                    │ PRAD tRF-QTL                │
                                    │ PRAD iso-eQTL               │
                                    │ PRAD sQTL                   │
                                    └──────────────┬──────────────┘
                                                   │
                                                   ▼
                                    ┌─────────────────────────────┐
                                    │      Data Validation        │
                                    └──────────────┬──────────────┘
                                                   │
                                                   ▼
                                    ┌─────────────────────────────┐
                                    │ Schema Harmonization        │
                                    │ SNP / Allele / Coordinate   │
                                    │ Genome Build                │
                                    └──────────────┬──────────────┘
                                                   │
                                                   ▼
                                    ┌─────────────────────────────┐
                                    │ Statistical Filtering       │
                                    └──────────────┬──────────────┘
                                                   │
                                                   ▼
                                          ┌────────┴─────────┐
                                          │                  │
                                          ▼                  ▼
                                    ┌─────────────┐    ┌───────────────┐
                                    │ GWAS ×      │    │ GWAS ×        │
                                    │ tRF-QTL     │    │ iso-eQTL/sQTL │
                                    └──────┬──────┘    └───────┬───────┘
                                           │                   │
                                           └─────────┬─────────┘
                                                     ▼
                                          ┌──────────────────────┐
                                          │ Candidate Regulatory │
                                          │ Loci                 │
                                          └──────────┬───────────┘
                                                     ▼
                                          ┌──────────────────────┐
                                          │ LD-based Integration │
                                          └──────────┬───────────┘
                                                     ▼
                                          ┌──────────────────────┐
                                          │ Colocalization       │
                                          └──────────┬───────────┘
                                                     ▼
                                          ┌──────────────────────┐
                                          │ Functional Analysis  │
                                          └──────────┬───────────┘
                                                     ▼
                                          ┌──────────────────────┐
                                          │ Clinical Validation  │
                                          └──────────┬───────────┘
                                                     ▼
                                          ┌──────────────────────┐
                                          │ External Validation  │
                                          └──────────────────────┘
```
## Data Sources
The project is designed to integrate several publicly available genomic datasets.

<b>1. Prostate Cancer GWAS</b>
Genome-wide association summary statistics for prostate cancer will be used to identify susceptibility variants and risk loci.
Required fields include:
```text
    - SNP
    - Chromosome
    - Position
    - Effect allele
    - Other allele
    - Effect size / OR
    - Standard error
    - P-value
    - Effect allele frequency
```

<b>2. Cancer-tRFQTL</b>
Cancer-tRFQTL data provide genetic associations between variants and tRNA-derived fragment expression across cancer types.
For this project, the prostate adenocarcinoma (PRAD) dataset will be prioritized.
Required information includes:
```text
    - SNP
    - Genomic position
    - Alleles
    - tRF identifier
    - Beta
    - P-value
    - FDR
    - MAF
    - Imputation quality
```

<b>3. Moradi et al. Prostate Cancer QTL Dataset</b>
The prostate cancer QTL dataset from Moradi et al. provides transcript-level regulatory associations.
```text
    The initial analysis will focus on:
    - cis-iso-eQTL
    - prostate-cancer-associated cis-QTL
    - trans-iso-eQTL
    - prostate-cancer-associated trans-QTL
```
The supplementary tables are expected to include:

```text   
    - SNP
    - Genomic position
    - Alleles
    - Isoform / splicing event
    - Beta
    - P-value
    - FDR
    - MAF
```

Splicing events may include:

```text 
    - cassette exon events
    - intron retention events
    - other alternative transcript events
```

<b>4. TCGA-PRAD</b>
TCGA prostate adenocarcinoma data may be used for downstream clinical and expression validation.
Potential analyses include:
```text 
    - expression differences,
    - Gleason score association,
    - tumor stage association,
    - survival analysis,
    - molecular subtype association.
```

<b>5. MINTbase / tRF Expression Resources</b>
tRF expression data may be used to validate whether prioritized tRFs are expressed in prostate cancer tissue.

<b>6. GTEx</b>
GTEx prostate tissue data may be used as an independent validation resource for:
```text
    - eQTL associations,
    - sQTL associations,
    - gene-level regulatory effects.
```

<b>7. LD Reference Data</b>

LD reference populations such as 1000 Genomes may be used for:
```text
    - proxy SNP identification,
    - LD expansion,
    - locus definition,
    - colocalization preparation.
```

<b> Project Structure </b>
```text
    pcatRFQTL/
    │
    ├── pyproject.toml
    ├── README.md
    ├── LICENSE
    ├── .gitignore
    ├── .env.example
    │
    ├── configs/
    │   ├── config.yaml
    │   ├── datasets.yaml
    │   └── thresholds.yaml
    │
    ├── data/
    │   ├── raw/
    │   │   ├── gwas/
    │   │   ├── trfqtl/
    │   │   ├── iso_eqtl/
    │   │   ├── sqtl/
    │   │   ├── tcga/
    │   │   ├── mintbase/
    │   │   ├── gtex/
    │   │   └── ld/
    │   │
    │   ├── interim/
    │   │   ├── harmonized/
    │   │   ├── filtered/
    │   │   ├── lifted/
    │   │   └── merged/
    │   │
    │   └── processed/
    │       ├── gwas/
    │       ├── trfqtl/
    │       ├── qtl/
    │       └── candidates/
    │
    ├── results/
    │   ├── tables/
    │   ├── statistics/
    │   ├── colocalization/
    │   ├── enrichment/
    │   └── validation/
    │
    ├── figures/
    │   ├── exploratory/
    │   ├── main/
    │   └── supplementary/
    │
    ├── notebooks/
    │   ├── 01_data_exploration.ipynb
    │   ├── 02_gwas_exploration.ipynb
    │   ├── 03_trfqtl_exploration.ipynb
    │   └── 04_candidate_exploration.ipynb
    │
    ├── scripts/
    │   ├── download_data.py
    │   ├── run_pipeline.py
    │   └── validate_data.py
    │
    ├── src/
    │   └── pcatrfqtl/
    │       ├── __init__.py
    │       ├── cli.py
    │       │
    │       ├── config/
    │       │   └── settings.py
    │       │
    │       ├── io/
    │       │   ├── readers.py
    │       │   ├── writers.py
    │       │   └── formats.py
    │       │
    │       ├── schemas/
    │       │   ├── gwas.py
    │       │   ├── trfqtl.py
    │       │   ├── eqtl.py
    │       │   └── sqtl.py
    │       │
    │       ├── preprocessing/
    │       │   ├── cleaning.py
    │       │   ├── harmonization.py
    │       │   ├── genome.py
    │       │   └── filtering.py
    │       │
    │       ├── gwas/
    │       │   ├── loader.py
    │       │   ├── filtering.py
    │       │   └── loci.py
    │       │
    │       ├── trfqtl/
    │       │   ├── loader.py
    │       │   ├── filtering.py
    │       │   └── annotation.py
    │       │
    │       ├── qtl/
    │       │   ├── iso_eqtl.py
    │       │   ├── sqtl.py
    │       │   └── trans_qtl.py
    │       │
    │       ├── integration/
    │       │   ├── overlap.py
    │       │   ├── ld.py
    │       │   ├── harmonize.py
    │       │   └── candidates.py
    │       │
    │       ├── colocalization/
    │       │   ├── preparation.py
    │       │   ├── coloc.py
    │       │   └── susie.py
    │       │
    │       ├── functional/
    │       │   ├── targets.py
    │       │   ├── enrichment.py
    │       │   └── pathways.py
    │       │
    │       ├── clinical/
    │       │   ├── tcga.py
    │       │   ├── survival.py
    │       │   └── associations.py
    │       │
    │       ├── validation/
    │       │   ├── gtex.py
    │       │   └── replication.py
    │       │
    │       ├── visualization/
    │       │   ├── plots.py
    │       │   ├── locus.py
    │       │   ├── volcano.py
    │       │   └── network.py
    │       │
    │       └── utils/
    │           ├── logging.py
    │           ├── validation.py
    │           └── statistics.py
    │
    └── tests/
        ├── test_gwas.py
        ├── test_trfqtl.py
        ├── test_harmonization.py
        ├── test_overlap.py
        └── test_statistics.py
```

