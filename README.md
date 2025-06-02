# Transcriptional Noise Analysis in Early *Drosophila* Embryo

This repository contains the code used to analyze transcriptional noise across different phenotypes as presented in the study:  
**"Mitotic reactivation and transcriptional bursting govern transcriptional noise in the early *Drosophila* embryo"**

## 📌 Overview

The analysis includes:
- Computing **mean**, **variance**, and **Fano Factor (FF)** of transcriptional activity
- Aggregating these metrics across **time windows** and **individual nuclei**
- Comparing noise statistics across multiple **phenotypes**

## 🚀 Pipeline Usage

To run the full analysis pipeline, use the following command in your Python environment:

```python
noise(input_path, ["phenotype_1/", "phenotype_2/"], output_folder)
