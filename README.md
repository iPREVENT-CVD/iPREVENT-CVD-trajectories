# Cardiovascular Risk Factors and Biomarkers in the Decades Prior to a Cardiovascular Event

## Table of Contents
- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [Authors](#authors)
- [Acknowledgements](#acknowledgements)
- [Citation](#citation)
- [Licence](#licence)

## About
- In this cohort study using individual-participant data pooled from four cohort studies, we evaluated sex- and age-specific trajectories of cardiovascular risk factors and circulating biomarkers in relation to cardiovascular events and over the adult life course.
- Related publications: https://doi.org/10.1093/eurjpc/zwaf764; https://doi.org/10.1186/s12916-023-02921-8; https://doi.org/10.1161/circulationaha.123.064386

## Installation
Step by step instructions to get your code running:
1. Clone the repository
git clone https://github.com/iPREVENT-CVD/iPREVENT-CVD-trajectories.git

2. Install R
Download and install R from https://www.r-project.org/
Recommended version: R 4.0+

3. Install RStudio (Recommended)
Download from https://www.rstudio.com/

4. Install required packages. At the start of the script required packages are stated. If these have not already been downloaded follow these steps to install packages prior to running scripts: (1) Open R or RStudio (2) run install.packages(c("package1", "package2"))

## Usage
The R code provided include 1) the linear and generalized linear mixed-effects models used to estimate average trajectories of cardiovascular risk factors and circulating biomarkers across individuals who experienced the primary outcome, individuals who survived and did not experience the primary outcome, and individuals who died from non-cardiovascular causes. analysis; and 2) the linear and generalized linear mixed-effects models used to evaluate sex-specific trajectories of cardiovascular risk factors and circulating biomarkers over the adult life-course. Additionally, two functions are provided that were used in the R scripts to run the models. For each analysis, an R script is available to create the figures as presented in the manuscript.

## Project Structure
```
├── function_timetoevent/              # R scripts
├── function_lifecourse/               # R scripts
├── 001_model_timetoevent/             # R scripts
├── 002_figure_creation_timetoevent/   # R scripts
├── 003_model_lifecourse/              # R scripts
├── 004_figure_creation_lifecourse/    # R scripts

└── README.md    # This file
```

## Dependencies
- R 4.0+
- RStudio (recommended)
- Packages required:
   - tidyverse
   - nlme
   - GLMMadaptive     

## Authors R code
- Marie de Bakker - University of Edinburgh
- Dave Yeung - University of Edinburgh
- Dorien Kimenai - University of Edinburgh
- A full list of the authors of the manuscript can be found here: placeholder

## Acknowledgements
This work was funded by the British Heart Foundation Intermediate Basic Science Research Fellowship (FS/IBSRF/23/25161) held by DMK, and MdB is supported by a Research Excellence Award from the British Heart Foundation (RE/24/130012). PW is supported by grant funding from the Chief Scientist Office, Scottish Government (HIPS/24/31). TMHE was supported by an Established Investigator E-Dekker grant (#03-002-2023-0036) and FIT-HEART consortium grant (#01-001-2024-0621) of the Dutch Heart Foundation. This Cardiovascular Health Study research was supported by NHLBI contracts HHSN268201200036C, HHSN268200800007C, HHSN268201800001C, N01HC55222, N01HC85079, N01HC85080, N01HC85081, N01HC85082, N01HC85083, N01HC85086, 75N92021D00006; and NHLBI grants U01HL080295, R01HL087652, R01HL103612, R01HL105756, R01HL120393, U01HL130114, and R01HL172803 with additional contribution from the National Institute of Neurological Disorders and Stroke (NINDS). Additional support was provided through R01AG023629 from the National Institute on Aging (NIA). A full list of principal CHS investigators and institutions can be found at CHS-NHLBI.org. The Trøndelag Health Study (HUNT) is a collaboration between HUNT Research Centre (Faculty of Medicine and Health Sciences, Norwegian University of Science and Technology NTNU), Trøndelag County Council, Central Norway Regional Health Authority, and the Norwegian Institute of Public Health. This Multi-Ethnic Study of Atherosclerosis research was supported by contracts 75N92025D00022, 75N92020D00001, HHSN268201500003I, N01-HC-95159, 75N92025D00026, 75N92020D00005, N01-HC-95160, 75N92020D00002, N01-HC-95161, 75N92025D00024, 75N92020D00003, N01-HC-95162, 75N92025D00027, 75N92020D00006, N01-HC-95163, 75N92025D00025, 75N92020D00004, N01-HC-95164, 75N92025D00028, 75N92020D00007, N01-HC-95165, N01-HC-95166, N01-HC-95167, N01-HC-95168 and N01-HC-95169 from the National Heart, Lung, and Blood Institute, and by grants UL1-TR-000040, UL1-TR-001079, and UL1-TR-001420 from the National Center for Advancing Translational Sciences (NCATS). Some mortality data were provided by the Bureau of Vital Statistics, New York City Department of Health and Mental Hygiene. The authors thank the other investigators, the staff, and the participants of the MESA study for their valuable contributions. A full list of participating MESA investigators and institutions can be found at https://mesa-nhlbi.org. This paper has been reviewed and approved by the MESA Publications and Presentations Committee. Support for the BNP dataset was provided by Roche Applied Sciences (Indianapolis, Indiana, USA). Funding support for the MESA family history dataset was provided by grant R21 HL81175. Funding support for the Renal Function dataset was provided by grant DK083538-01. This Whitehall II Study research was supported by a British Heart Foundation Programme grant (RG/16/11/32334). We thank all of the participating civil service departments and their welfare, personnel, and establishment officers, the British Occupational Health and Safety Agency, the British Council of Civil Service Unions, all participating civil servants in the Whitehall II study, and all members of the Whitehall II Study team at UCL and Oxford. MdB and DMK had full access to all the data in the study and takes responsibility for the integrity of the data and the accuracy of the data analysis.

## Citation
If you use this code in your research please cite:
placeholder

## Licence
This project is licensed under the [MIT Licence](LICENSE)

## Contact
For questions or issues please contact:
- Dorien Kimenai, email: dorien.kimenai@ed.ac.uk

