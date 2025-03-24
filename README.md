![Snakemake](https://img.shields.io/badge/snakemake-Workflow-blue.svg?style=flat-square) 

[![Cookiecutter](https://img.shields.io/badge/built_with-Cookiecutter-green.svg?style=flat-square)](https://cookiecutter.readthedocs.io/)



# Pipeline for Enrichment Analysis of Variants in Case Cohorts

This pipeline evaluates the enrichment of variants with **adjustable parameters**, e.g.:
- **Allele Frequency (AF) < 0.01**
- **CADD scores > 20** and / or
- **spliceAI scores > 0.2**

in a case cohort compared to the gnomAD reference cohort. Tresholds for which to filter for can be adapted. It is particularly beneficial for small case cohorts where statistical power is insufficient to detect variants in disease-causing genes. Analyzing large individual-level variant datasets (e.g., UK Biobank) would ideally improve control cohort size, but requires substantial computational resources and data storage infrastructure. Additionally, existing datasets mostly represent Caucasian populations, limiting genetic studies of understudied populations. This can be circumvented by this approach. 

### **Advantages**
- **Enhanced Statistical Power**: Integrates ~68,000 genome samples from the gnomAD database (caucasian genomes) as a control cohort, significantly increasing the number of control samples.
- **Broad Variant Analysis**: Enables comparison of rare, potentially deleterious variants across case and reference populations.

### **Limitations**
- **Lack of Individual-Level Data**: Since individual sample-level data from gnomAD genomes are unavailable, this information is lost, potentially inflating differences due to multiple variants in the same individual.
- **No Linkage Disequilibrium Accounting**: The pipeline cannot account for linkage disequilibrium between variants.



---

### **How to Run**
1. Prepare a list of gene names or a bed file with regions you wish to analyze
2. Create a conda environment based on the yaml file provided and get containers for **VEP**, **BCFTOOLS** and **BEDTOOLS** from repositories.
3. Retrieve reference files necessary to run the pipeline (see config).
4. Adjust values in the config to your analysis.
5. Run snakemake workflow in the following way: *snakemake --profile scripts/snakemake_profile_cookiecutter*

The workflow is written to run in a slurm based job submisson environment, but can also be run locally by altering the profile CMD line argument in the snakemake call. Values for the partition to submit to and arguments to bind certain folders might have to be adjusted / added by the user.


### **Acknowledgments**
gnomAD

snakemake

cookiecutter

