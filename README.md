# Description
This repository contains data documentation and code supporting the analyses in "Post-Acute Outcomes After Early Gabapentinoid–Opioid Use in Skilled Nursing Facilities".

## Repository Contents
- `data_documentation/`: Contains general study information and additional supplementary material.
- `code/`: Example data management and analysis programs.
- `LICENSE`: The license under which this repository is shared.
- `README.md`: This file, providing an overview of the repository.

## Data Documentation
The `data_documentation/` directory contains the following files:
- `Data_Documentation.xlsx`: Includes the following tabs:
  -  `OVERVIEW`: Overview of Excel tabs.
  -  `Project Overview`: Overview of authors, study question, and data sources used.
  -  `Neuropathic Pain`: ICD-10-CM codes for the neuropathic pain covariate.
  -  `Spine Pain & Degen Disord`: ICD-10-CM codes used for the spine pain/degenerative disorders covariate.
  -  `Chronic Pain Other`: ICD-10-CM codes used for the chronic pain due to other causes covariate.
  -  `Other Musculo Pain`: ICD-10-CM codes used for the other musculoskeletal pain covariate.

## Code
The `code/` directory contains the following programs:
- Key data preparation programs (`00_data_prep/`):
  - `00_create_artificial_censoring_vars.sas`: Create variables for dates of censoring due to treatment strategy deviation.
  - `01_prep_follow_up_wide_ds.sas`: Prepare follow-up variables for outcome-specific wide-format datasets.
  - `02_prep_person_period_ds.sas`: Create outcome-specific person-period datasets for pooled logistic regression models.
- Key analytic programs (`01_analysis/`):
  - `00_artificial_censoring_weight_estimation.R`: Estimate inverse probability of remaining uncensored weights for treatment strategy deviation using baseline covariates.
  - `01_outcome_estimation.R`: Fit weighted pooled logistic regression models to estimate outcome risks, risk differences, and risk ratios.
  - `02_bootstrap_estimates.R`: Generate and save estimates from bootstrap replicates for subsequent calculation of 95% credible intervals.
  - `03_construct_ci.R`: Combine point estimates with 95% credible intervals derived from bootstrap estimates and tabulate results at 30 and 99 days after SNF admission.

These programs illustrate key data preparation and analysis steps rather than a complete, executable pipeline. Additional code can be shared upon reasonable request.

## Funding
Research reported in this GitHub repository was supported by the National Institute on Aging of the National Institutes of Health under award number R01AG077620. This project was supported by project funding totaling $3,524,470, of which 100% is financed with federal funds and 0% is financed by nongovernmental sources. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
