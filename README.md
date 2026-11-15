# PGScaler

Explore the [PGS Catalog](https://www.pgscatalog.org): search the polygenic scores (PGS / PRS)
of one or more EFO / MONDO traits, gather every published evaluation of these scores and
produce a report that helps choose one.

## What it does

1. **Search scores**: every score associated with each trait id in the PGS Catalog (and its
   child traits), or a given score id (`PGSxxxxxx`). Traits are grouped: child traits with
   their queried parent, and any explicit group from the config. The rest of the pipeline
   works per group (one report section per group). Scores outside `date_min`–`date_max` or
   with fewer than `n_snp` variants are dropped at this step (their evaluations and scoring
   files are not downloaded); the report header lists the filters and the groups left empty.
2. **Collect evaluations**: every performance metric reported for these scores (one PGS
   Catalog *PPM* = one evaluation of one score, in one sample, for one phenotype).
3. **Classify each evaluation**: is it about the queried trait? Was it published
   independently from the score? Are the metrics for the score alone or for a model with
   covariates? (see [How evaluations are classified](#how-evaluations-are-classified))
4. **Check key loci** (optional): download the scoring files and look for given variants
   (e.g. APOE for Alzheimer's disease).
5. **Summarize** each score in a decision table, give it an evidence level, list the studies
   that compared several scores head to head, draw figures and write an HTML report.

## Installation

The environment is managed with [pixi](https://pixi.sh) (dependencies in
[pyproject.toml](pyproject.toml): pandaspgs, hydra, matplotlib, plotly…):

```bash
cd PGScaler
pixi install
```

## Usage

```bash
# from the parent directory of the repo (any working directory works)
pixi run --manifest-path PGScaler/pyproject.toml python PGScaler/src/pgscaler/main.py

# override any option with Hydra
pixi run --manifest-path PGScaler/pyproject.toml python PGScaler/src/pgscaler/main.py \
    'trait=[EFO_0004340,MONDO_0005090]' target_ancestry=EAS

# change or disable the filters
pixi run --manifest-path PGScaler/pyproject.toml python PGScaler/src/pgscaler/main.py \
    date_min=2020 date_max=null n_snp=1000
```

Outputs always go to `outputpath`, resolved from the repository root when relative
(`outputs/` by default), whatever the directory the command is run from. Hydra's own run
directories and log files are disabled in the config.

API answers are cached by pandaspgs, so a second run is fast. Scoring files are only
downloaded for key-loci checks, once per score (results cached in `outputs/cache/`).
If no score passes the filters, the run stops with a message and no report is written.

## Configuration

[src/config/config.yaml](src/config/config.yaml):

| Option | Description |
|---|---|
| `outputpath` | Output directory; a relative path is relative to the repository root (default `outputs`). |
| `trait` | List of EFO / MONDO / HP trait ids, or PGS ids (`PGS002332`). Find ids on the PGS Catalog or [EBI OLS](https://www.ebi.ac.uk/ols4). |
| `include_child_traits` | Also fetch the scores of the child traits of each trait, as the PGS Catalog website does (default `true`). E.g. *bipolar disorder* also brings the scores of *bipolar I disorder* and *bipolar II disorder*. The `source_traits` column says under which trait(s) each score was found. |
| `group_child_traits` | Merge a queried trait into another queried trait it descends from in the ontology (default `true`): *bipolar I disorder* + *bipolar II disorder* + *bipolar disorder* -> one *bipolar disorder* section. The top-most queried ancestor wins. |
| `trait_groups` | Explicit groups `{name: [trait ids]}`, for traits that are not parent/child in the ontology or to choose the name: `depression: [MONDO_0002050, MONDO_0002009]`. Override the automatic grouping (a trait listed alone in a group stays separate). |
| `target_ancestry` | Ancestry of the population the score will be applied to: `EUR`, `EAS`, `SAS`, `AFR`, `AMR`, `GME`, `ASN`. Used for the `target_evaluated` column, notes and sorting. |
| `date_min`, `date_max` | Filter: keep only the scores published between these years (inclusive). `null` disables a bound. |
| `n_snp` | Filter: keep only the scores with at least `n_snp` variants (removes scores with too few variants). `null` disables it. Scores with an unknown publication date or number of variants are dropped when the corresponding filter is set. |
| `eval_patterns` | Optional, per trait id or group name (patterns of the traits of a group are combined): case-insensitive regex that recognizes the evaluations of the trait itself among the evaluated phenotypes. Use it when automatic matching is wrong. |
| `key_loci` | Optional, per trait id or group name: variants to look for in the scoring files (GRCh38 positions). |

Example:

```yaml
outputpath: "outputs"

trait:
  - MONDO_0004985   # bipolar disorder (its child traits are fetched and grouped with it)
  - EFO_0009963     # bipolar I disorder -> merged into bipolar disorder
  - MONDO_0002050   # depressive disorder  \ grouped as "depression"
  - MONDO_0002009   # major depressive disorder /
  - MONDO_0004975   # Alzheimer disease
  - PGS002332       # a single score

include_child_traits: true
group_child_traits: true
trait_groups:
  depression: [MONDO_0002050, MONDO_0002009]

target_ancestry: EUR

date_min: 2019      # filters, null to disable
date_max: 2026
n_snp: 100

eval_patterns:
  MONDO_0004975: "^(?!.*biomarker).*(alzheimer|\\bAD\\b)"  # AD, excluding biomarker studies
  MONDO_0005090: "schizo|psychos"                          # also accept psychosis

key_loci:
  MONDO_0004975:
    - name: APOE
      chr: 19
      snps: {rs429358: 44908684, rs7412: 44908822}   # rsid: GRCh38 position
      window: 500000                                  # +/- bp counted as "region"
```

## Outputs

| File | Content |
|---|---|
| `report.html` | Self-contained HTML report ([structure](#report)). |
| `scores_summary.csv` | One row per (trait group, score) with every decision criterion ([columns](#decision-table-scores_summarycsv)). |
| `performance_metrics.csv` | Long table: one row per (evaluation, metric), with sample, ancestry, cohorts, covariates and the `primary` / `external` flags. Start here to go beyond the summary. |
| `head_to_head.csv` | Scores evaluated side by side in the same study, sample and phenotype, ranked. |
| `figures/` | PNG figures of each trait group, also embedded in the report ([list](#figures)). |
| `cache/key_loci.csv` | Key-loci check results. |

### Report

`report.html` is a single file (≈ 6 MB, plotly.js and figures embedded): it opens offline
and can be shared as is. Light theme, in English.

1. **Header**: date, target ancestry, links to each trait group, filters applied and the
   groups left without any score by the filters.
2. **All scores**: interactive plotly figure of every kept score, publication date × trait
   group, colored by evidence level, marker size growing with the number of variants (log
   scale, rescaled between the smallest and largest score shown). Hovering a point shows the PGS id, name, reported trait,
   publication, number of variants, level, evaluations, ancestries and caveats; clicking it
   opens the score page on the PGS Catalog.
3. **How to read this report**: short guide (evidence levels, comparability of metrics).
4. **One section per trait group**, with its trait ids (linked to the PGS Catalog):
   - key figures: scores, evaluations of the trait, externally validated scores, scores
     evaluated in the target ancestry, level A scores;
   - shortlist: table of the level A and B scores;
   - head-to-head comparisons (collapsible, one per study, sample and metric);
   - figures;
   - all scores of the group (collapsible table).

### Figures

Per trait group, in `figures/` (file names start with the group key):

| File | Content |
|---|---|
| `<group>_overview.png` | Left: number of variants (log) vs. number of evaluations on the trait, one point per score, size ∝ √(largest number of cases, or N/10), colored by evidence level. Right: scores per publication year, stacked by evidence level. The left panel is left out when no score of the group has been evaluated on the trait. |
| `<group>_OR.png`, `<group>_HR.png` | Forest plots: one point per primary evaluation with its confidence interval, log scale. Colored by evaluation ancestry (EUR, EAS, AFR, other / mixed); filled = external evaluation, hollow = by the authors. |
| `<group>_PGSAUROC(nocovariates).png` | Same, AUC of the score alone (reference line at 0.5). |
| `<group>_AUROC_C-index.png` | Same, AUC / C-index of models with covariates. |
| `<group>_PGSR2(nocovariates).png` | Same, R² of the score alone. |

A forest plot is only drawn when at least one evaluation reports the metric.

## How evaluations are classified

**Primary evaluation** (`primary`): the evaluated phenotype is the queried trait. Scores are
often evaluated on other phenotypes (a body-weight PGS tested on hearing difficulty, an
Alzheimer PGS on amyloid burden or cognition). These evaluations are kept in
`performance_metrics.csv` but **excluded from every metric of the summary**; they are only
counted (`n_evals_other_traits`) and listed (`other_phenotypes`).

- With `eval_patterns`: the phenotype matches the regex.
- Without: the phenotype shares a word stem (first 5 letters, generic words such as
  *disease*, *disorder*, *risk* ignored) with the trait label or the score's reported trait,
  or contains the acronym of the trait label (*major depressive disorder* -> *MDD*).

Automatic matching is a heuristic: check `other_phenotypes` and add a pattern if needed.
Close phenotypes (dementia vs. Alzheimer's disease, MCI, age at onset) are **not** matched
unless the pattern says so.

**External evaluation** (`external`): published by another publication than the one that
developed the score. It is the best protection against overfitting and optimistic
reporting, but it may still come from the same consortium or reuse the same cohorts.

**Score alone vs. full model**: the PGS Catalog stores metrics as reported by the authors.
An AUROC or R² is considered to be for the score alone only when the authors say so
(*PGS AUROC (no covariates)*, *PGS R2 (no covariates)*) or when no covariate is reported;
otherwise it is a full-model metric.

## Metrics

### Metrics reported in the PGS Catalog

| Metric | Type of trait | What it measures | How to read it |
|---|---|---|---|
| **AUROC** (AUC) | Binary | Discrimination: probability that a random case has a higher predicted risk than a random control. | 0.5 = random, 1 = perfect. Usually for **score + covariates**: age and sex alone often give 0.6–0.8, so a high AUC says little about the score itself. |
| **PGS AUROC (no covariates)** | Binary | AUC of the score alone. | The cleanest measure of how well the score discriminates. PGS rarely exceed 0.6–0.7. |
| **Incremental AUROC** | Binary | AUC gain when the score is added to a covariate model. | What the score brings beyond known factors. Often 0.01–0.05. |
| **C-index** | Time to event | AUC equivalent for survival data (who develops the disease first). | Same reading as AUROC. |
| **OR** (odds ratio) | Binary | Odds of disease associated with the score, from a logistic regression. | 1 = no effect. **Check the unit**: per standard deviation (most common; PGS typically 1.2–2), per allele, or top vs. bottom quantile (much larger values). Not comparable across units. |
| **HR** (hazard ratio) | Time to event | Same as OR for an incident outcome (Cox model). | Same reading and unit caveats as OR. |
| **β** | Continuous (sometimes binary) | Regression coefficient of the score. | Depends on the scale of the phenotype and of the score: only comparable within a study. |
| **R²** | Continuous (pseudo-R² for binary) | Phenotype variance explained, usually by **score + covariates**. | 0–1. Same caveat as AUROC: covariates can explain most of it. |
| **PGS R2 (no covariates)** | Continuous | Variance explained by the score alone. | The cleanest measure for quantitative traits. |
| **Incremental R²** | Continuous | R² gain when the score is added to a covariate model. | What the score brings beyond covariates. |
| Other metrics | | NRI, Δ C-index, AUPRC, OR/HR between extreme quantiles, cumulative incidence… | Kept in `performance_metrics.csv`, not summarized. |

Metrics are as reported by the authors: the PGS Catalog does not recompute them. Always
check covariates (`covariates` column), the unit of the effect and the exact phenotype
before comparing two values.

### Decision table (`scores_summary.csv`)

Every metric column is computed over the **primary evaluations only**, across all their
samples and ancestries (an AUC measured in EUR and one in AFR count the same).

| Column | Definition |
|---|---|
| `query_trait`, `trait_label` | Trait group key (root trait id or explicit group name) and its label. |
| `pgs_id`, `name`, `trait_reported` | Score id, name given by its authors, trait as reported by the authors. |
| `source_traits` | Trait id(s) the score was found under (queried trait or one of its child traits). |
| `method`, `n_variants` | Construction method (P+T, LDpred2, PRS-CS…) and number of variants. Scores of a few dozen variants work with a genotyping array; genome-wide scores (10⁵–10⁶ variants) need good imputation and coverage. |
| `dev_publication`, `dev_date`, `dev_journal`, `dev_doi` | Publication that developed the score (first author and year), its date, journal and DOI. |
| `gwas_ancestry`, `dev_ancestry` | Ancestry distribution of the GWAS the weights come from, and of the sample used to tune the score. A score usually performs worse in other ancestries. |
| `n_evals` | Number of primary evaluations. |
| `n_evals_other_traits`, `other_phenotypes` | Evaluations on other phenotypes (excluded from metrics) and the first 4 of them. |
| `n_external_pubs`, `external_by` | Number of independent publications with a primary evaluation, and their first authors. |
| `eval_ancestry`, `target_evaluated` | Ancestries of the primary evaluations; whether `target_ancestry` is among them. |
| `max_eval_n`, `max_eval_cases` | Size of the largest evaluation sample (total and cases; 0 cases = continuous trait or not reported). |
| `auc_pgs_only` | **Max** AUROC of the score alone. |
| `auc_with_covariates` | **Max** AUROC / C-index of a model with covariates. |
| `or_median`, `hr_median` | **Median** OR / HR over evaluations (median rather than max: mixed units make extremes meaningless). |
| `r2_pgs_only`, `r2_with_covariates` | **Max** R² of the score alone / of the full model. |
| `<locus>_snps`, `<locus>_region_variants` | Key-loci check: key SNPs found in the score (`none` if absent) and number of score variants within the window. A locus can be absent because it is modelled separately (e.g. APOE genotype as a covariate): check the publication. |
| `harmonized_grch38`, `scoring_file` | A harmonized GRCh38 scoring file exists; link to the original scoring file. |
| `evidence` | Evidence level (below). |
| `notes` | Automatic warnings: never evaluated externally, never in the target ancestry, AUC only with covariates, < 200 cases, evaluations on other phenotypes, no harmonized file. |

Max values are optimistic summaries (best sample, best model): use them to spot a score,
then read the underlying evaluations in `performance_metrics.csv` or the report.

### Evidence levels

Based on primary evaluations only. "Large" = ≥ 1,000 cases for a binary trait, or N ≥ 10,000
when no cases are reported (continuous trait).

| Level | Rule | Meaning |
|---|---|---|
| **A** | At least one independent publication **and** a large evaluation | Replicated, reasonably precise. |
| **B** | Independent publication **or** large evaluation | Promising, but either possibly optimistic (only by its authors) or imprecise. |
| **C** | Only small evaluations by the developers | Weak evidence. |
| **D** | No evaluation of the trait itself | Only tested on related phenotypes, or not at all. |

Within a trait, scores are sorted by evidence level, then evaluated in the target ancestry,
then number of independent publications, then size of the largest evaluation. The level
does not rank performance: two A scores must be compared on their metrics.

### Head-to-head comparisons

Metrics from different studies are **not comparable**: samples, case definitions, ancestry,
covariates and units differ. The only fair comparison is between scores evaluated in the
same publication, phenotype, ancestry, cohorts and metric. `head_to_head.csv` and the report
list these groups and rank the scores in each (rank 1 = highest value; for every metric
summarized here, higher means better). A benchmark such as Monti et al. 2024, which compares
seven methods built on the same GWAS in three independent biobanks, is often the most useful
piece of information for choosing.

## How to choose a score

1. Keep scores of the right trait (check `trait_reported`, not only the queried id).
2. Prefer level A, then B; read `notes`.
3. Check `target_evaluated` and `gwas_ancestry` against your population.
4. Compare candidates with head-to-head comparisons, otherwise with metrics of the score
   alone (`auc_pgs_only`, `r2_pgs_only`), never with full-model metrics across studies.
5. Check feasibility: `n_variants` vs. your genotyping/imputation, `harmonized_grch38`,
   key loci, and the license and terms of use on the score page.
6. Read the publication of the chosen score: unit of the effects, covariates, overlap
   between your cohort and the samples used to build the score (e.g. UK Biobank).

## Code structure

| Module | Role |
|---|---|
| [main.py](src/pgscaler/main.py) | Entry point (Hydra): runs the pipeline and writes the outputs. |
| [catalog.py](src/pgscaler/catalog.py) | PGS Catalog access: traits and their children, trait grouping, scores, filters, evaluations, key-loci checks in the scoring files. |
| [analysis.py](src/pgscaler/analysis.py) | Classification of evaluations, decision table, evidence levels, head-to-head comparisons. |
| [plots.py](src/pgscaler/plots.py) | Interactive plotly figure and matplotlib figures. |
| [report.py](src/pgscaler/report.py) | HTML report. |
| [config.yaml](src/config/config.yaml) | Configuration. |

## Notes

- `pandaspgs.get_scores(trait_id=…)` does not return the scores of child traits: they are
  fetched separately when `include_child_traits` is set.
- `pandaspgs.get_performances()` crashes when no sample of an evaluation lists a cohort:
  evaluations are read from the raw REST answers instead (same endpoint and cache).
- Rerunning with new traits only queries what is not cached yet.
