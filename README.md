# Microbiota Developmental Trajectory Analysis for Preterm Infants

## Overview

This repository contains the code and data analysis for a study investigating the role of gut microbiota maturation in preterm infants and its relationship with late-onset sepsis (LOS). Our findings suggest that antibiotic exposure disrupts microbiota development, leading to a deficiency in bacterial enzymes such as DL-endopeptidase, which are important for the activation of the host’s immune system. Specifically, we identify a deficiency in DL-endopeptidase-producing *Enterococcus faecium* as a key factor in delayed microbiota maturation and increased LOS risk. The analysis tools in this repository allow for the exploration of microbiota development, NOD2 activation, and the effect of antibiotics on microbiota composition in preterm infants.

## Code Structure

This repository includes the following key analysis scripts:

- **data_tidy**: Data preprocessing and cleaning functions to prepare microbiome data for analysis.
- **profile**: Functions for profiling microbial communities using 16S rRNA amplicon sequencing.
- **maaslin2_diff**: Differential abundance analysis using MaAsLin2 to identify microbiota features associated with clinical outcomes.
- **shap_value**: Calculation of Shapley values to assess feature importance in predictive models.
- **cox_hr_sur_curve**: Survival analysis using Cox proportional hazards models to explore the impact of microbiota features on late-onset sepsis risk.
- **NOD2_analysis**: Analysis of NOD2 activation and its relationship to microbiota composition.
- **cluster_trend**: Cluster analysis to identify trends in microbiota maturation across cohorts.
- **slow_fast_sepsis_ab**: Analysis of microbiota composition in relation to slow vs. fast microbiota maturation and sepsis risk.
- **mediation_analysis**: Mediation analysis exploring the causal relationship between antibiotic exposure, microbiota maturity, and sepsis incidence.
- **Poilt_metagenome**: data analysis pipeline for the clinical intervention.
- **wet_valid**: Display and analysis of animal experiment data to validate findings from human cohorts.

## Purpose of the Analysis

The goal of these analyses is to:

- **Characterize gut microbiota development** in preterm infants across multiple cohorts from China, the USA, and the UK.
- **Identify biomarkers** (such as DL-endopeptidase) associated with delayed microbiota maturation and increased risk of late-onset sepsis.
- **Provide a predictive model** to assess microbiota maturation trajectory, which could serve as a tool to predict LOS risk and guide clinical interventions.
- **Assess the impact of early-life antibiotic exposure** on microbiota maturation and its consequences on immune function and disease susceptibility.

## How to Use the Code

Clone this repository to your local machine:
   ```bash
   git clone https://github.com/ZJJY-Bioinformatics/infant_gut_MM.git
   ```


