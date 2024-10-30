# MDR1 Variants

This repository contains the code used to run the analyses described in the paper "Computational Analysis of MDR1 Variants Predicts Effect on Cancer Cells via their Effect on mRNA Folding" by Tal Gutman and Tamir Tuller.

## Overview

MDR1 encodes for p-glycoprotein (p-gp), a transmembrane efflux pump that plays a crucial role in patient response to chemotherapy. This project investigates three key variants in the MDR1 gene (T1236C, T2677G, and T3435C) that are suspected to impact p-gp function.

Through multiple computational analyses, we examine how these variants affect different steps of MDR1 expression to better understand:
- Their mechanism of action
- Impact on chemotherapy response 
- Effect on patient survival

## Repository Structure

- `MDR1_code/`: Contains the main analysis notebooks and scripts
- `co_trans_code/`: Contains the co-translational folding model implementation (unpublished)

## Installation

This project requires three separate conda environments due to different dependencies:

1. Main environment:
```bash
conda create --name mdr1_main --file MDR1_requirements.txt
```

2. Splicing analysis environment:
```bash 
conda create --name mdr1_splice --file spliceAI_requirements.txt
```

3. Enformer analysis environment:
```bash
conda create --name mdr1_enformer --file Enformer_requirements.txt
```















