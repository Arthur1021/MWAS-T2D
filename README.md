## Analysis workflow

### 1. meQTL analysis and population-specific meta-analysis

**Directory:**

```text
01_meQTL_and_meta/
```

This step identifies genetic variants associated with DNA methylation levels and combines meQTL association results across available datasets within each population.

**Main scripts:**

- `process.with.pca.general.sh`  
  Performs meQTL association analysis with population-specific principal components and other covariates.

- `generate.metal.script.pl`  
  Generates METAL scripts for population-specific meQTL meta-analysis.

- `generate.metal.script.for.singles.pl`  
  Generates METAL scripts for CpGs or datasets requiring separate processing.

Population-specific CpG-SNP association estimates are combined using METAL with:

```text
SCHEME STDERR
```

This corresponds to fixed-effect inverse-variance weighted meta-analysis using effect estimates and standard errors.

Typical output includes:

```text
CpG
SNP
effect allele
other allele
beta
standard error
P value
```

These results are used for DNA methylation prediction model development.

---

### 2. DNA methylation prediction model training, selection, and validation

**Directory:**

```text
02_train_models_selection_validation/
```

This step develops population-specific genetic prediction models for DNA methylation.

**Model training**

```text
1_train.R
```

Trains genetically predicted methylation models using cis-meQTL information.

**Prepare validation genotype data**

```text
2_prepare_selection_validation_genotype.R
```

Prepares genotype data for model selection and validation.

**Population-specific model selection**

```text
3_selection_validation_AA_parallel.R
3_selection_validation_CA_parallel.R
3_selection_validation_EA_parallel.R
3_selection_validation_HA_parallel.R
```

Performs parallelized model selection and validation in each population.

**Generate final model files**

```text
4_generate_model_files.R
```

Generates the final population-specific methylation prediction model files used for downstream MWAS.

**Prediction and validation**

```text
5_1_generate_predicted_methylation.R
5_2_generate_PCA_for_each_ethnic.R
5_3_generate_covariate.R
5_4_validation_rsq.R
```

These scripts generate genetically predicted methylation values, calculate population-specific principal components, prepare covariates, and calculate validation R².

---

### 3. MWAS association testing and cross-population meta-analysis

**Directory:**

```text
03_association_test_and_meta/
```

**Population-specific association testing**

```text
01_association.sh
```

Performs population-specific methylome-wide association testing of genetically predicted DNA methylation with type 2 diabetes.

**Cross-population meta-analysis**

```text
02_meta.R
```

Combines population-specific CpG-T2D association statistics across populations.

The cross-population MWAS meta-analysis is based on population-specific Z statistics and sample size.

**Bias and inflation correction**

```text
03_Controlling_bias_and_inflation.R
```

Applies the `bacon` method to evaluate and correct bias and inflation in MWAS test statistics.

Outputs include both original and bacon-corrected association statistics.

---

### 4. Conditional analysis

**Directory:**

```text
04_conditional_test/
```

Conditional analyses evaluate whether MWAS signals remain associated with T2D after conditioning the T2D GWAS summary statistics on known T2D index variants.

A common set of **1,289 independent T2D index variants** from the multi-ancestry T2D GWAS was used for conditioning.

Population-matched 1000 Genomes reference panels were used for LD estimation:

```text
AA  -> AFR
AsA -> EAS
EA  -> EUR
HA  -> AMR
```

**Prepare conditional GWAS statistics**

```text
03_prepare_conditional_GWAS.py
```

Prepares population-specific conditioned T2D GWAS summary statistics.

**Re-run MWAS using conditional GWAS**

```text
04_MWAS_using_conditional_GWAS.R
```

Repeats MWAS using the conditioned GWAS results.

**Extract conditionally significant signals**

```text
05_extract_sig_after_conditional_MWAS.R
```

Identifies CpGs that remain significant after conditioning.

**Extract signals no longer significant after conditioning**

```text
06_extract_non_sig_after_conditional_MWAS.R
```

**Evaluate non-European-only findings at known T2D loci**

```text
07_nonEuropean_only_established_T2D_loci.R
```

Quantifies the empirical contribution of the multi-population design.

Non-European-only CpGs are defined as CpGs reaching Bonferroni significance in AA, AsA, and/or HA but not in EA.

These CpGs are mapped to the **611 known T2D loci** reported in the multi-ancestry T2D GWAS.

In the current analysis:

- **399** CpG-T2D associations were detected only through non-European population analyses.
- **367 of 399 (91.98%)** mapped to **82 unique known T2D loci**.

---

### 5. Mendelian randomization analysis

**Directory:**

```text
05_MR_analysis/
```

MR analysis is used as a supportive sensitivity analysis for MWAS-significant CpGs.

Because the MR analyses use genetic instruments derived from the same underlying meQTL resources and T2D GWAS used for the MWAS, they should not be interpreted as fully independent replication.

**Main scripts:**

```text
01_MR_analysis.R
02_summary.R
```

For each MWAS-significant CpG:

1. methylation-associated variants are selected;
2. ancestry-matched LD clumping is performed;
3. variants are harmonized with the corresponding population-specific T2D GWAS.

Population-matched 1000 Genomes LD reference panels are used.

Depending on the number of available instruments:

- Wald ratio is used for a single instrument.
- Inverse-variance weighted MR is used for multiple instruments.

Additional sensitivity analyses are performed where applicable.

---

### 6. Colocalization analysis

**Directory:**

```text
06_coloc_analysis/
```

Colocalization analysis evaluates whether the methylation association and the T2D GWAS association within a region are consistent with a shared underlying genetic signal.

**Extract cis regions**

```text
01_extract_sig_cis_region.R
```

Extracts regions surrounding MWAS-significant CpGs.

**Prepare GWAS data**

```text
02_prepare_GWAS.py
```

Prepares population-specific T2D GWAS variants within the corresponding regions.

**Run colocalization**

```text
03_run_coloc.R
coloc_pipeline.R
```

Colocalization is performed using the `coloc` R package.

Posterior probabilities are calculated for the standard five hypotheses:

```text
H0: neither trait is associated
H1: only methylation is associated
H2: only T2D is associated
H3: both are associated but with different causal variants
H4: both are associated and share a causal variant
```

The primary statistic is:

```text
PP.H4.abf
```

A value of:

```text
PP.H4.abf > 0.70
```

is considered supportive evidence of colocalization.

---

### 7. Transcriptome-wide association analysis

**Directory:**

```text
07_TWAS/
```

This directory contains or is reserved for transcriptome-wide association analysis used to evaluate whether genes linked to MWAS findings also show genetically predicted expression associations with T2D.

Some TWAS analyses depend on external reference models and controlled-access data.

---

### 8. Proteome-wide association analysis

**Directory:**

```text
08_PWAS/
```

PWAS is used to evaluate protein-level evidence for genes prioritized through the methylation analysis.

**Main scripts:**

```text
01_get_gene_list.R
02_extract.R
```

These scripts obtain genes corresponding to prioritized MWAS findings and extract the corresponding PWAS results.

---

### 9. GO and KEGG enrichment analysis

**Directory:**

```text
09_GO_KEGG/
```

**Main script:**

```text
code.R
```

Performs functional enrichment analysis of prioritized genes using:

- Gene Ontology (GO)
- Kyoto Encyclopedia of Genes and Genomes (KEGG)

---



