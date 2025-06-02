# Transcriptional Noise Analysis in Early *Drosophila* Embryo

This repository contains the code used to analyze transcriptional noise across different phenotypes as presented in the study:  
**"Mitotic reactivation and transcriptional bursting govern transcriptional noise in the early *Drosophila* embryo"**

## Overview

The analysis includes:
- Computing **mean**, **variance**, and **Fano Factor (FF)** of transcriptional activity
- Aggregating these metrics across **time windows** and **individual nuclei**
- Comparing noise statistics across multiple **phenotypes**

## Pipeline Usage

To run the full analysis pipeline, use the following command in your Python environment:

```python
noise(input_path, ["phenotype_1/", "phenotype_2/"], output_folder)

## Arguments:
input_path: Root path where phenotype folders are stored

["phenotype_1/", "phenotype_2/"]: List of phenotype folder names (with trailing slashes)

output_folder: Path where you want the output figures and results to be stored

## Folder Structure (per phenotype):
Each phenotype folder must contain:
