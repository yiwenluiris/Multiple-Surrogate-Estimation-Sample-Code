# Leveraging Error-Prone Algorithm-Derived Phenotypes: Enhancing Association Studies for Risk Factors in EHR Data

This is the code for our paper "Leveraging Error-Prone Algorithm-Derived Phenotypes: Enhancing Association Studies for Risk Factors in EHR Data"

## Overview
![Figure_abstract](https://github.com/yiwenluiris/Multiple-Surrogate-Estimation-Sample-Code/assets/93690224/603a4f0e-2e57-43cf-ae12-03b9b26ab2d5)



## Description
**Objectives:** It has become increasingly common for multiple computable phenotypes from electronic health records (EHR) to be developed for a given phenotype. However, EHR-based association studies often focus on a single phenotype. In this paper, we develop a method aiming to simultaneously make use of multiple EHR-derived phenotypes for reduction of bias due to phenotyping error and improved efficiency of phenotype/exposure associations.<br />

**Materials and Methods:** The proposed method combines multiple algorithm-derived phenotypes with a small set of validated outcomes to reduce bias and improve estimation accuracy and efficiency. The performance of our method was evaluated through simulation studies and real-world application to an analysis of colon cancer recurrence using EHR data from Kaiser Permanente Washington.<br />

**Results:** In settings where there was no single surrogate performing uniformly better than all others in terms of both sensitivity and specificity, our method achieved substantial bias reduction compared to using a single algorithm-derived phenotype. Our method also led to higher estimation efficiency by up to 30% compared to an estimator that used only one algorithm-derived phenotype. <br />

**Discussion:** Simulation studies and application to real-world data demonstrated the effectiveness of our method in integrating multiple phenotypes, thereby enhancing bias reduction, statistical accuracy and efficiency.<br />

**Conclusions:** Our method combines information across multiple surrogates using a statistically efficient seemingly unrelated regression framework. Our method provides a robust alternative to single-surrogate-based bias correction, especially in contexts lacking information on which surrogate is superior.

## Code and sample data
**Code** Run _multi_surrogate_code.R_ with sample data or generic data in the same format to try the method<br />
**Sample data** The generated synthetic data _sample_data.csv_ is provided in the repo. 

## Questions?
Feel free to contact _yiwenlu at sas.upenn.edu_ for any questions regarding this repo. 
