# EEG Signal Analysis for Alzheimer’s Biomarker Detection using Phase Locking Value during Olfactory Stimulation

This repository contains the code, report, and supplementary materials for our project on analyzing EEG signals to identify biomarkers for Alzheimer's Disease (AD) using Phase Locking Value (PLV) during olfactory stimulation. This project was completed as part of the Signals and Systems course under the supervision of **Prof. Hamid Aghajan** at Sharif University of Technology.

## Project Overview

### Objective
The objective of this project is to investigate the use of Phase Locking Value (PLV) as a biomarker for Alzheimer's Disease (AD) by analyzing EEG signals during olfactory stimulation. The study explores the difference in PLV between AD patients and healthy individuals when exposed to specific odors, using EEG data.

### Methodology
- **Data Collection**: EEG data was collected from participants during olfactory stimulation using lemon (frequent) and rose (rare) odors.
- **Data Preprocessing**: The EEG signals were preprocessed using a standard pipeline, including filtering, artifact removal, and re-referencing, to ensure clean and usable data.
- **Phase Locking Value (PLV) Calculation**: PLV was calculated between specific EEG channels (Fz and Cz) to assess the synchronization of neural activity during the olfactory task.
- **Statistical Analysis**: Box plots, Gaussian fitting, and p-value calculations were performed to determine the significance of PLV differences between AD patients and healthy controls.

## Files in this Repository

- `Project_Report.pdf`: The comprehensive report detailing the project's objectives, methods, and findings.
- `SignalsSystems_Project1402.pdf`: The project documentation including task definitions, data descriptions, and methodology.
- `NikooMoradi_400101934.m`: MATLAB script for preprocessing EEG signals and calculating PLV.
- `PreProcessing_maincode.m`: MATLAB script for the complete EEG preprocessing pipeline.
- `codes`: Folder containing all the scripts for pipeline.
- `Dataset`: Folder containing different parts of the dataset.

## Authors
- **Nikoo Moradi**

## Supervisor
- **Prof. Hamid Aghajan**, Department of Electrical Engineering, Sharif University of Technology
