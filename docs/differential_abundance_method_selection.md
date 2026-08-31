# Differential-abundance analysis: training decision framework

## Current limitation

The present dataset contains three simulated samples: two labelled T2D and one labelled Control. It is suitable for testing file structures and transformations, but it cannot support differential-abundance inference, estimation of group variance, reliable covariate adjustment, or false-discovery-rate control.

## Provisional strategy for a future T2D–microbiome study

### Primary analysis: ANCOM-BC2

Use ANCOM-BC2 on the untransformed integer count table with a prespecified model and filtering strategy. ANCOM-BC2 provides bias correction for compositional microbiome data and supports covariates, continuous exposures, multiple groups and repeated-measures designs. Its reported abundance effects should be interpreted as model-based estimates under the method's assumptions, not as direct measurements of absolute microbial load.

Do not provide ANCOM-BC2 with percentages or manually calculated CLR values. Preserve raw counts and let the method perform its own estimation and bias correction.

### Sensitivity analysis

Choose at least one method with a different modelling strategy:

- **MaAsLin2** when the principal need is flexible multivariable association modelling, including mixed effects and complex metadata;
- **ALDEx2** for a compositional Monte Carlo analysis with standardized effect sizes, particularly for a simpler comparison with adequate biological replication;
- **LinDA** as a scalable CLR-based, bias-corrected linear-model sensitivity analysis, including designs requiring covariates or mixed effects.

Agreement across methods strengthens robustness but does not prove biological truth. Disagreement should be reported and investigated rather than resolved by selecting the method producing the most significant taxa.

## Method comparison

| Method | Main input | Principal strength | Important limitation or caution | Proposed role |
|---|---|---|---|---|
| ANCOM-BC2 | Raw integer counts and metadata | Bias correction; covariates; multi-group and repeated-measures analyses | Results depend on model assumptions, filtering, zero handling and adequate replication | Primary analysis |
| MaAsLin2 | Counts or transformed abundance table, depending on the specified workflow | Flexible multivariable and mixed-effects modelling | Normalization and transformation choices affect the estimand and results | Multivariable sensitivity analysis |
| ALDEx2 | Raw counts and design information | Dirichlet Monte Carlo sampling; log-ratio framework; effect sizes | Can be conservative or underpowered; reference/denominator and replication matter | Compositional sensitivity analysis |
| LinDA | Count table with CLR-based modelling | Fast, flexible linear models with compositional bias correction | Zero replacement and model assumptions require scrutiny | Alternative sensitivity analysis |
| DESeq2 | Raw counts | Mature negative-binomial framework for RNA-seq | Not designed specifically for microbiome compositionality; benchmark performance is dataset-dependent | Not the default primary microbiome method |
| Naive tests on relative percentages | Percentages | Simple | Ignore closure, sparsity, confounding and multiplicity | Do not use as primary inference |

## Prespecified estimand and model

For a cross-sectional T2D study, define the primary estimand before analysis, for example:

> Adjusted difference in the abundance of each prespecified genus between adults with T2D and controls in the target population.

The exact interpretation will depend on the selected method. Report effect estimates, uncertainty intervals where available, raw p-values and multiplicity-adjusted q-values. Do not report only lists of significant taxa.

## Covariate strategy for T2D studies

Potential covariates include age, sex, recruitment site, sequencing batch, recent antibiotics, stool consistency, diet and relevant medications. Metformin deserves particular attention because treatment can be strongly related to T2D status and microbiome composition.

Covariates must be selected using the scientific question and a causal framework. Do not adjust automatically for every measured variable. BMI, diet, medication and metabolic biomarkers may be confounders, mediators or consequences depending on the estimand. A causal diagram and a prespecified primary model should guide adjustment, followed by sensitivity analyses.

## Required safeguards

1. Complete contamination assessment and technical QC before differential-abundance testing.
2. Define the taxonomic level, prevalence filter and library-size criteria before examining group p-values.
3. Verify sample independence or model repeated measures explicitly.
4. Examine influential samples, batch effects and group imbalance.
5. Control multiplicity across tested features and contrasts.
6. Report null findings and method disagreement.
7. Validate important signals in an independent dataset or targeted assay where feasible.
8. If absolute abundance is scientifically important, incorporate microbial-load information such as spike-ins, flow cytometry or quantitative PCR rather than inferring it from relative sequencing alone.

## Training decision

No differential-abundance test will be performed on the current three-sample table. The next practical step is to verify the R environment, then use a clearly labelled larger simulated dataset to learn ANCOM-BC2 and a sensitivity method. Results from that dataset will remain pedagogical and will not be interpreted as evidence about T2D.

## Primary references and documentation

- Lin H, Peddada SD. ANCOM-BC2: [Nature Methods article](https://doi.org/10.1038/s41592-023-02092-7) and [official Bioconductor manual](https://bioconductor.org/packages/release/bioc/manuals/ANCOMBC/man/ANCOMBC.pdf).
- Mallick H et al. MaAsLin2: [PLOS Computational Biology](https://doi.org/10.1371/journal.pcbi.1009442).
- ALDEx2: [official Bioconductor vignette](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html).
- Zhou H et al. LinDA: [Genome Biology](https://doi.org/10.1186/s13059-022-02655-5).
- Nearing JT et al. Cross-dataset comparison of differential-abundance methods: [Nature Communications](https://doi.org/10.1038/s41467-022-28034-z).
