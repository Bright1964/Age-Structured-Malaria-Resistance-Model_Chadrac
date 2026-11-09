# Age-structured malaria model: ITN control and insecticide resistance in the DRC

**Code and simulated-data workflow for the manuscript**

*An age-structured model for assessing the impact of vector control on malaria transmission and insecticide resistance evolution*

This repository is prepared from the authors' R scripts to organise the analyses underlying **five main figures and Supplementary Figure S1**. All code comments and figure outputs are in English. The repository is **not yet a verified, complete reproducibility release**: see [Important limitations](#important-limitations-before-public-release).

## Figure-to-code map

| Figure in manuscript | Analysis and panels | Generation code | Exported plot data |
|:--|:--|:--|:--|
| **Figure 1** | A/B: resistance invasion index; C/D: malaria basic reproduction numbers at mosquito equilibria, for two net products | `scripts/01_simulate_main.R` and `scripts/02_figures_main.R` | `output/Figure_1_invasion_data.csv`, `output/Figure_1_R0_data.csv` |
| **Figure 2** | Provinces (A: 9, B: 9, C: 8), 2009–2024, projected resistance and campaign markers | `scripts/04_figure2_provinces.R` | `output/Figure_2_monthly_provincial_resistance.csv`, `output/Figure_2_campaign_markers.csv` |
| **Figure 3** | A: prevalence, B: resistance allele frequency; C: proportions of cases/deaths averted, comparing products at 50% coverage | `scripts/01_simulate_main.R` and `scripts/02_figures_main.R` | `output/Figure_3_monthly_data.csv`, `output/Figure_3_impact_data.csv` |
| **Figure 4** | A/B: pyrethroid-only; C/D: next-generation; prevalence and resistance at 25%, 50%, 75% coverage | Same two main scripts | `output/Figure_4_monthly_data.csv` |
| **Figure 5** | A/B: proportions of malaria cases and deaths averted at 25%, 50%, 75% coverage | Same two main scripts | `output/Figure_5_impact_data.csv` |
| **Figure S1** | Age-group-specific cases and deaths averted, by net type and coverage | `scripts/03_figureS1.R` (uses main model) | `output/Figure_S1_plot_data.csv`, `output/age_specific_cases_deaths_averted.csv` |

Corresponding PDF/PNG files are created under `figures/`. The 26-province figure uses the author-supplied Excel workbook at `data/bednet_distribution_by_province_2008_2022.xlsx`; a long-form CSV and campaign-events CSV are also included. **Rights to publish these programme records on a public repository must be confirmed.**

## Reproduction workflow

1. Install a recent R version and RStudio (optional). Install the R dependencies by running `source("install_packages.R")` from the project root.
2. The original provincial campaign workbook is now included at `data/bednet_distribution_by_province_2008_2022.xlsx`. Check `data/README.md` for provenance, file hash, data validation, and public-sharing permissions.
3. Run the project from the repository root:

   ```r
   source("run_all.R")
   ```

4. Inspect the six figure files in `figures/` and the exact input tables for each figure in `output/`.
5. Record the software versions in `output/R_sessionInfo.txt` and compare the generated plots **numerically and visually** with the submitted manuscript before releasing a public archive.

### Fast trial (not for publication)

`FAST_MODE=TRUE` reduces the number of intervention draws and the spin-up/calibration iterations; its figures are not substitutes for final publication runs.

```r
Sys.setenv(FAST_MODE = "TRUE")
source("run_all.R")
```

For publication, restart R, then run:

```r
Sys.setenv(FAST_MODE = "FALSE")
source("run_all.R")
```

The **daily numerical step is fixed at 1 day**. The code uses one-day shifts along the human and mosquito infection-age dimensions. Increasing `DT` to 2 or 5 without redesigning those shifts would no longer solve the same numerical scheme. For the threshold curves of Figure 1 only, coverage is sampled every **0.01 (1 percentage point)** rather than 0.005 (0.5 percentage point) to shorten the sweep. It does not modify the underlying model equations.

The supplied provincial file is now included and passes structural checks (26 provinces, 2009–2024, all 416 campaign flags 0/1, 103 positive flags). If you remove it, `run_all.R` will skip Figure 2. **These checks do not verify the historical accuracy of the programme records, nor prove the 2019 PBO assumption.**

## Files

```text
.
├── run_all.R
├── install_packages.R
├── scripts/
│   ├── 01_simulate_main.R       # Age-structured daily dynamic simulations; threshold sweep
│   ├── 02_figures_main.R        # Figures 1, 3, 4, 5 only
│   ├── 03_figureS1.R            # Age-specific simulation and Figure S1
│   └── 04_figure2_provinces.R   # One provincial mechanistic model; Figure 2 only
├── data/                       # Original Figure 2 workbook, derived campaign CSVs, and data notes
├── docs/
│   ├── figure_mapping.csv
│   └── reproducibility_notes.md
├── figures/                    # Generated plots (ignored by Git until released)
├── output/                     # Generated numerical files (ignored by Git until released)
└── cache/                      # Local computational cache (never release by default)
```

## Methods and interpretation

The human compartments are stratified by seven age groups; infectious age is tracked by day. Mosquitoes are stratified by `SS`, `RS`, and `RR` genotypes. Simulations cover **12 years** for the prospective interventions and **2009–2024** for the provincial historical campaign scenarios. The prospective scenarios compare pyrethroid-only and generic next-generation ITNs at 25%, 50% and 75% peak coverage, with three-year campaigns. All scenarios start at a resistant-allele frequency of 0.05.

The revised plotting code uses **raw model-derived monthly summaries**, removing the older display-only smoothing and forced end-frequency trends. The Figure 3 and Figure 5 bar centres use **medians** (instead of means) to match the manuscript captions, with 2.5th and 97.5th percentiles across the parameter draws. Pre-intervention prevalence is based on the actual calibrated endemic state rather than an imposed 50% plotted baseline. Consequently, these new plots may **not exactly match the current manuscript**; the revised manuscript and/or numerical assumptions must be checked before submission.

## Important limitations before public release

- **Figure 2 input now included:** the workbook contains 26 provinces and 2009–2024 mass-distribution indicators, with no missing or invalid binary values. Its filename refers to 2008–2022 despite its 2009–2024 columns. The provincial `PBO_FIRST_YEAR=2019` rule still classifies later campaigns as PBO **by assumption**, because the original Excel file records campaign occurrence but not product type. These projections are not observed resistance allele frequencies. Source provenance and redistribution rights have not been independently verified.
- **Uncertainty provenance:** some intervention parameter draws in `make_draws()` are sampled from assumed truncated normal distributions, not read from a public posterior file. Describe these as **scenario/parameter uncertainty intervals**, unless calibrated posterior inputs and their provenance are supplied.
- **Unvalidated numerical workflow:** R itself was unavailable in the preparation environment, so the refactored scripts have been reviewed structurally but **not run to completion**. Numerical reproducibility, timings, and fit to all manuscript panels must be verified in RStudio.
- **Simulation and publication agreement:** changing from display-only smoothed trajectories to raw model trajectories, and plotting medians rather than means, may change figure appearance/numeric centres. Check all results and revise figure descriptions as required; never manually adjust numerical predictions to match the submitted graphics.
- **Methodological assumptions:** `BURDEN_CALIBRATION_FACTOR=0.30`, the simplified efficacy of generic next-generation nets, the use of `pmax(0, no-ITN minus ITN)` for prevented burden, and parameter bounds should all be explained or justified in the manuscript. A no-ITN arm that does worse is truncated at zero by the original code.
- **Version/licensing:** obtain coauthor agreement on the code/data licence and verify the rights to redistribute the original programme records. This package intentionally does not assign a licence or invent a permanent repository URL.

Review `docs/reproducibility_notes.md` for a checklist before sending the revision to the journal.

## Data/code availability wording (insert only after a real public release)

> The R scripts and simulated numerical datasets supporting Figures 1–5 and Supplementary Figure S1 are available at **[INSERT VERIFIED PUBLIC REPOSITORY URL/DOI]**. The repository contains figure-specific reproducibility scripts, simulation settings, source data tables and execution instructions. Any restrictions on underlying provincial ITN programme data are described in the repository.

Do not cite this placeholder as if it were a published repository. Consider freezing the verified GitHub commit in **Zenodo** and citing its assigned DOI in the manuscript.

## Static repository check

The script `docs/validate_structure.py` can check code delimiter balance, manuscript figure mapping and the absence of the old redundant plot definitions. This is **not** a substitute for a full run in R.

```bash
python docs/validate_structure.py
```
