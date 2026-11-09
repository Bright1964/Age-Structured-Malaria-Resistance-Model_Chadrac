# Provincial ITN mass-distribution campaigns (Figure 2)

## Supplied input

`bednet_distribution_by_province_2008_2022.xlsx` is an **unmodified copy** of the Excel workbook provided by the authors for Figure 2. Its file name ends in `2008_2022`, but its actual year columns run from **2009 through 2024**.

The workbook has one sheet (`Sheet1`), with **26 distinct provinces** in its first column (`province`), followed by **16 annual fields** (`2009` through `2024`). All 416 province-year cells contain binary campaign indicators (1 = campaign recorded, 0 = no campaign recorded); none are missing. There are **103 positive indicators**, meaning 103 recorded province-year campaign occurrences, not 103 distinct national campaigns.

`province_campaigns_2009_2024.csv` is a long-form, machine-readable export of all 416 cells. `campaign_events_2009_2024.csv` contains just the 103 province-year pairs flagged as `1`. The data in these two CSVs were transcribed programmatically from the attached workbook **without changing values**. The project's R script reads the original Excel file, not these derived CSVs.

SHA-256 checksum of the supplied Excel file:

`a446f5e5b93c3046246d1f6cf492ec54317c4d364351469108047cdee0aa00ab`

## How the data enter Figure 2

From the repository root, `source("run_all.R")` invokes `scripts/04_figure2_provinces.R`, which reads the original workbook using `readxl::read_excel()`. Province row order is preserved and the 26 provinces are divided into three figure panels (9/9/8). The output CSVs are:

- `output/Figure_2_monthly_provincial_resistance.csv`: predicted **resistance allele frequency**, not observed insecticide-resistance measurements.
- `output/Figure_2_campaign_markers.csv`: event years and assigned campaign products used for the plot.

The file contains **only mass-distribution indicators**, not distribution volumes, actual net coverage, nor net product type. The mechanistic model imposes a **50% peak coverage scenario** at each campaign; this coverage is not inferred from the workbook. It labels campaigns from 2019 onward as **PBO (assumed)** and uses generic next-generation ITN effectiveness as a proxy: this is an **unverified scenario assumption**, not a recorded provincial product history.

## Provenance and release permission

These data were provided for this research-code preparation, but no original source document, NMCP provenance metadata or public-redistribution licence accompanied the workbook. The historical entries were checked for **format and internal consistency only**, not independently verified against original campaign records. **Confirm that the data owner permits publishing the workbook and its derived CSVs before pushing the repository to a public GitHub account.** If permission is withheld, do not publish these files; describe the restrictions honestly in the manuscript, and discuss reviewer access with the editor.
