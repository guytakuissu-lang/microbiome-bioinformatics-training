# Prespecified ANCOM-BC2 repeated-simulation training protocol

## Status and purpose

This document prespecifies the next ANCOM-BC2 training module before the larger simulated datasets are generated or analyzed.

The purpose is to test the reproducibility of the computational workflow and to learn how filtering, covariate adjustment, sparsity, and repeated simulation affect apparent signal recovery. It is not a formal validation or benchmark of ANCOM-BC2. No simulated result will be interpreted as biological evidence about type 2 diabetes (T2D), metformin, or the gut microbiome.

## Primary training question

Under a prespecified multivariable simulation with more than 50 retained taxa, how consistently does the existing ANCOM-BC2 workflow recover known artificial T2D group effects while controlling discoveries among simulated group-null taxa?

## Secondary training questions

1. How often are true artificial group signals removed by prevalence filtering?
2. How different are end-to-end signal recovery and recovery conditional on a taxon being retained?
3. How biased and variable are estimated group log effects across simulation replicates?
4. Does adjustment separate an artificial metformin effect from the conditional T2D group contrast?
5. Are conclusions sensitive to the prespecified prevalence threshold?

## Study design

### Base-case sample size

- Total samples: 120
- Control: 60
- T2D: 60

This sample size is a training choice. It is not justified as a power calculation for a real study.

### Number of taxa

- Simulated taxa: 100
- Target retained after the primary filter: at least 80
- Minimum acceptable retained in any fitted replicate: more than 50

If a replicate retains 50 or fewer taxa, it will be flagged as failing the base-case QC target. It will not be silently discarded. The cause will be investigated and the replicate will remain represented in the QC report.

### Pilot replication

- Number of simulation replicates: 20
- Replicate identifiers: `R001` to `R020`
- Base seed: `20260901`
- Replicate seed: base seed plus the replicate number

Twenty replicates are sufficient to test pipeline automation but insufficient for stable method-performance claims. Any estimated operating characteristics will be labelled pilot estimates with substantial Monte Carlo uncertainty.

## Simulated metadata

Each replicate will contain the following variables:

| Variable | Prespecified distribution or allocation | Reference level |
|---|---|---|
| `group` | 60 Control and 60 T2D | Control |
| `age` | Continuous adult age, generated within a prespecified plausible range | Centred mean age |
| `sex` | Balanced Female/Male allocation within group | Female |
| `metformin` | No Control exposure; 30 untreated and 30 treated participants within T2D | No |
| `batch` | Balanced Batch1/Batch2 allocation within group | Batch1 |
| `library_size` | Variable positive integer sequencing depth | Not applicable |

Age generation will use the same distribution in both groups in the base case so that age imbalance is not introduced unintentionally.

### Important metformin limitation

Metformin exposure occurs only in the T2D group. The design matrix remains estimable because treated and untreated participants both occur within T2D, but there are no metformin-exposed Control participants. The adjusted `groupT2D` coefficient is therefore a conditional contrast anchored to metformin = No. It is not a marginal T2D effect and does not solve the causal positivity limitation that would arise in a real observational study.

## Count-generation model

The simulation will generate raw integer counts. Percentages and CLR values will not be used as ANCOM-BC2 input.

For each replicate:

1. assign each taxon a fixed baseline log abundance;
2. include a mixture of common, intermediate, and rare taxa;
3. add the prespecified taxon-specific covariate effects on the log scale;
4. add sample-level taxon variability to avoid deterministic compositions;
5. transform latent log abundances to sample-specific multinomial probabilities;
6. draw the observed count vector using the simulated library size.

This produces compositional counts by construction. Zeros may arise naturally for low-probability taxa. The base case will not add a separate structural-zero or zero-inflation mechanism. Such mechanisms may be considered only in a separately prespecified future scenario.

## Library-size generation

Library sizes will vary across samples and remain positive integers. The generator will be calibrated to produce a realistic training range without creating a systematic library-size difference between T2D and Control groups.

Each replicate must report:

- minimum, median, and maximum library size;
- group-specific library-size summaries;
- any sample with zero total counts;
- agreement between generated library sizes and count-table column totals.

## Prespecified artificial effects

Effect-bearing taxa will be assigned before simulation and held constant across the 20 replicates.

### T2D group effects

- Total artificial group signals: 10 taxa
- Positive effects: 5 taxa
- Negative effects: 5 taxa
- Target absolute log-effect magnitude: approximately 0.6 in the base case
- Common/intermediate group signals: 8 taxa
- Rare group signals: 2 taxa

Including two rare group signals allows the workflow to distinguish failure of statistical detection from removal during filtering.

### Other covariate effects

Artificial effects will also be assigned to non-overlapping taxa for:

- age per 10-year increase;
- Male versus Female;
- metformin Yes versus No;
- Batch2 versus Batch1.

These effects will be disjoint from the 10 group-effect taxa in the base case. This simplifies interpretation of group recovery and prevents an artificial taxon from carrying both a group and medication effect.

### Simulated group-null taxa

All taxa without an assigned T2D group effect are group-null for performance calculations, even if they carry an age, sex, metformin, or batch effect.

## Primary filtering rule

- Detection definition: count greater than or equal to 10
- Minimum prevalence: detected in at least 20% of samples
- Library-size filter inside ANCOM-BC2: none beyond prior QC in this training module

The primary filter is applied without reference to group p-values or q-values.

Each replicate will retain a taxon-level filter report containing:

- total count;
- detection prevalence;
- retained/excluded status;
- simulated effect status for every covariate.

## Sensitivity filtering rules

After the primary analysis has been completed, the same simulated replicates may be reanalyzed with minimum prevalence thresholds of 5% and 50%, while retaining the count >= 10 detection definition.

These sensitivity thresholds are prespecified here. They will not replace the 20% primary threshold based on which produces more desirable discoveries.

## Primary analysis model

The same adjusted ANCOM-BC2 formula will be fitted in every eligible replicate:

```text
group + age_centered + sex + metformin + batch
```

Prespecified settings:

- raw integer count input;
- taxa in rows and samples in columns;
- BH multiplicity adjustment;
- pseudocount sensitivity enabled;
- structural-zero assessment by group enabled;
- one computational worker initially, to preserve deterministic training behaviour;
- alpha = 0.05.

The primary artificial T2D discovery definition is:

```text
diff_robust_groupT2D == TRUE
```

The corresponding effect estimate, standard error, p-value, BH-adjusted q-value, pseudocount-sensitivity status, and simulated truth will be retained for every analyzed taxon.

## Prespecified QC checks

Every replicate must pass or explicitly report the following checks:

1. count values are finite, non-negative integers;
2. sample identifiers are unique and identically ordered in counts and metadata;
3. taxon identifiers are unique and match the truth table;
4. count-table column totals equal the generated library sizes;
5. metadata variables have no unexpected missing values;
6. group, sex, batch, and treatment allocations match the prespecified design;
7. the adjusted design matrix is full rank;
8. more than 50 taxa remain after the primary filter;
9. all artificial effects are assigned to the intended taxa and covariates;
10. all output tables contain the expected ANCOM-BC2 coefficient columns;
11. failed, warned, or non-converged model fits are recorded rather than silently omitted.

## Primary performance summaries

Performance summaries apply to the `groupT2D` coefficient only unless explicitly stated otherwise.

### End-to-end true-positive proportion

The denominator contains all 10 simulated group-effect taxa, including those removed by filtering.

```text
end-to-end TPR = robustly detected true group signals / 10
```

This metric captures both filter loss and model non-detection.

### Analysis-conditional true-positive proportion

The denominator contains only true group-effect taxa retained and tested in the replicate.

```text
conditional TPR = robustly detected retained true signals / retained true signals
```

This metric must not be reported without the end-to-end metric because conditioning on retention can hide signal loss caused by filtering.

### False-discovery proportion

```text
FDP = robust discoveries among simulated group-null taxa / max(total robust discoveries, 1)
```

The mean FDP across replicates is an empirical pilot estimate related to false-discovery control. With only 20 replicates, it will not be presented as a definitive FDR estimate.

### False-positive proportion

```text
FPP = robust discoveries among retained simulated group-null taxa / retained group-null taxa
```

### Estimation accuracy

For retained true group-effect taxa:

- signed bias: estimated minus simulated log effect;
- mean absolute error;
- root mean squared error;
- sign agreement;
- approximate interval coverage when the required uncertainty output is available.

Approximate intervals based on estimate +/- 1.96 standard errors will be labelled approximate and will not be treated as validated confidence intervals for method benchmarking.

### Filtering loss

For each replicate:

- total taxa removed;
- true group signals removed;
- group-null taxa removed;
- prevalence distributions for retained and removed taxa.

## Across-replicate summaries

For each primary metric, report:

- all 20 replicate-level values;
- mean and median where appropriate;
- standard deviation or interquartile range;
- minimum and maximum;
- number of successful, warned, and failed model fits.

Plots must show replicate-level variation rather than only a single aggregate value.

## Interpretation rules

1. Do not call an undetected simulated signal absent.
2. Do not call a filtered signal statistically non-significant; it was not tested.
3. Do not interpret a null simulated effect as proof that a corresponding biological pathway is absent.
4. Do not select the prevalence threshold that maximizes desired discoveries.
5. Do not describe average FDP from 20 replicates as a definitive FDR result.
6. Do not generalize base-case findings to high sparsity, structural zeros, repeated measures, low biomass, contamination, or real T2D cohorts.
7. Do not interpret artificial metformin or group effects as biological findings.

## Planned outputs

The implementation should generate, at minimum:

```text
data/simulated_ancombc_large_design.tsv
data/simulated_ancombc_large_truth.tsv
results/ancombc2_repeated_simulation_qc.tsv
results/ancombc2_repeated_simulation_taxon_results.tsv
results/ancombc2_repeated_simulation_metrics.tsv
results/ancombc2_repeated_simulation_summary.tsv
figures/ancombc2_repeated_simulation_performance.png
```

Large replicate-specific count tables should not automatically be committed. A deterministic generator, seeds, compact design/truth files, and summarized outputs may provide better reproducibility with lower repository size. The final versioning policy will be decided after measuring generated file sizes.

## Implementation sequence

1. version this protocol;
2. implement one base-case dataset generator;
3. run QC and visually inspect the single large dataset;
4. run ANCOM-BC2 on the single dataset and verify coefficient extraction;
5. generalize the validated code to 20 replicates;
6. compute prespecified metrics;
7. generate figures and inspect them visually;
8. document deviations from this protocol;
9. only then consider additional scenarios or methods.

## Criteria for completing the pilot module

The pilot is complete when:

- all 20 seeds are reproducible;
- every replicate has a documented QC and model status;
- no failed replicate is silently removed;
- primary metrics are calculated exactly as prespecified;
- filter loss is separated from model non-detection;
- figures show replicate-level uncertainty;
- the README and methods notes clearly state that this is a training pilot, not a method benchmark.

## Future scenarios requiring separate prespecification

The following are deliberately outside the base case:

- stronger group or batch imbalance;
- metformin exposure in both clinical groups;
- overlapping biological effects on the same taxa;
- high zero inflation or structural zeros;
- repeated measures and random effects;
- low-biomass contamination;
- taxonomic misclassification;
- alternative abundance distributions;
- comparison with MaAsLin2, ALDEx2, or LinDA;
- absolute-abundance information from spike-ins, qPCR, or flow cytometry.

Each extension should have a separate rationale, estimand, simulation mechanism, and analysis plan.
