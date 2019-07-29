# Malnutrition

Collaboration to study the effect of malnutrition on pneumonia in LMICs

# Protocol

This protocol will be updated as decisions regarding the analysis are made.

## Definitions of malnutrition and categories of malnutrition severity

Malnutrition definitions

low weight for age - acute malnutrition
low weight for height - acute wasting
low height for age - stunting

Malnutritions severities

'>-2SD / z score or >80% - normal
<-2SD to <-3SD / z score or 60-80% - moderate
<-3SD / z score or <60% - severe

## Missing data

Study authors contacted where there was either missing data or data was not provided in a stratified fashion for malnutrition mortality.

Where only the odds ratio, but not the original counts from the 2x2 tables were provided, we estimated these from the confidence intervals as per (http://onlinelibrary.wiley.com/doi/10.1002/sim.2287/abstract). This was required for one study (Nathoo 1993).


## Modelling choices

Since the studies [frequently collapsed categories of malnutrition severity](</Data/Collapsed_Uncollapsed numbers.csv>) , we chose to estimate the odds ratios using a Bayesian model, as this will allow us to use all of the data to estimate the odds ratios for each possible comparison.

The modelling can handle the fact that, for some studies, the event data for the malnutrition categories are collapsed. We need to make a decision how to handle the fact that the total number in each malnutrition category is also collaped.

In [01_maln.md](Scripts/01_maln.md) I have explored some options for dealing with this problem.

## Key Points Regarding Specific Studies

Tupasi 1988 - This study uses an upper cut off of 75% of expected weight for age which is different to the upper cutoff of 80% used in most studies

Ballard 1995 - This study uses upper cut offs of 90% expected weight for age and weight for height and an upper cut off of 95% expected height for age which is different to the upper cutoff of 80% used in most studies

Sehgal - study previously included but excluded as includes bronchiolitis in cases

Nathoo - new study included as able to calculate original numbers from relative risk with CIs 

Tupasi, 1990 - two publications from same study so only the publication with most complete data used

Zabihullah and Webb - study authors indicate that many children who were severely ill did not have their height recorded. This is likely to introduce bias which underestimates low w/h as a risk factor for pneumonia mortality. This may be something relevant to other studies.

Lazzerini/Hooli/Enarson - these three studies cover the same children on the same database all looking at weight/age. Only Lazzerini included as more recent and includes more children. Hooli also had lots of missing data making its results slightly questionable.

Nantanda - definition of malnutrition not included in published study but clarified by study author via e-mail

Kuti - Malnutrition case fatality numbers calculated from original dataset provided by study author against WHO reference growth standards

## Uploaded files

- Incidence_estimates.csv contains the results of studies investigating incidence of pneumonia episodes in relation to malnutriton
- Mortality_Numbers.csv contains the number of children who survived or died in relation to malnutrition status in each study
- Mortality_odds_ratios.csv contains odds ratios calculated from data within studies in terms of mortality risk in relation to malnutrition status


## Preliminary results

In the scripts folder if you click on the file [01_maln.md](Scripts/01_maln.md) you will find that the results of the initial modelling.
