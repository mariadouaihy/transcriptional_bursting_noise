# Transcriptional Noise Analysis in Early *Drosophila* Embryo

This repository contains the code used in the study:

**"Mitotic reactivation and transcriptional bursting govern transcriptional noise in the early *Drosophila* embryo"**

---

## Goal

This code quantifies transcriptional noise across different phenotypes by calculating:

- **Mean**, **variance**, and **Fano Factor (FF)** of transcriptional activity
- Computed both:
  - On **sliding time windows**
  - **Per nucleus**
- Comparison of transcriptional noise between **multiple phenotypes**

---

## Folder Structure

To run the code, data for each phenotype must be stored in a standardized format:
input_path/
├── phenotype_1/
│ ├── rawData/ # Contains raw MS2 fluorescence time series
│ └── DeconvOutput/ # Output from the burst deconvolution model
├── phenotype_2/
│ ├── rawData/
│ └── DeconvOutput/


---

##  Preprocessing Requirement

Before running this pipeline, the raw signal must be **deconvolved** using the model from this repository:

👉 [https://github.com/mariadouaihy/stochastic_transcriptional_model](https://github.com/mariadouaihy/stochastic_transcriptional_model)

> This model estimates the transcriptional bursts from raw MS2 traces.  
> Its output should be saved in the `DeconvOutput/` subfolder.

---

##  How to Run

Use the main function:

```python
noise(input_path, ["phenotype_1/", "phenotype_2/"], output_folder)
