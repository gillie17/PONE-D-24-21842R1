*** UK Biobank -- colorectal cancer dataset creation ***

version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

use "CR_new_extract.dta", clear // updated dataset

label variable cr_cox_censor10_age "Censor 10 years for Cox analysis"
label variable cr_c_cox_censor10_age "Censor 10 years for Cox analysis"
label variable cr_r_cox_censor10_age "Censor 10 years for Cox analysis"


/* MERGE CLINICAL DATA */

merge 1:1 n_eid using "CR_clinical_1.dta"
drop if _merge==2 
drop _merge

merge 1:1 n_eid using "CR_clinical_2.dta"
drop if _merge==2 
drop _merge

merge 1:1 n_eid using "CR_clinical_3.dta"
drop if _merge==2 
drop _merge


/* MERGE PRS */

merge 1:1 n_eid using "CR_140PRS_risks.dta"
drop if _merge==2 
drop _merge


/* MERGE BLOOD BIOMARKERS */

merge 1:1 n_eid using ///
	"C:\Users\gillian.dite\Dropbox\GTG\UKB Data\UKB_Biomarkers.dta"
drop if _merge==2 
drop _merge

* drop the second instance fields
drop n_30600_1_0 n_30610_1_0 n_30620_1_0 n_30630_1_0 n_30640_1_0 n_30650_1_0 ///
	n_30660_1_0 n_30670_1_0 n_30680_1_0 n_30690_1_0 n_30700_1_0 n_30710_1_0 ///
	n_30720_1_0 n_30730_1_0 n_30740_1_0 n_30750_1_0 n_30760_1_0 n_30770_1_0 ///
	n_30780_1_0 n_30790_1_0 n_30800_1_0 n_30810_1_0 n_30820_1_0 n_30830_1_0 ///
	n_30840_1_0 n_30850_1_0 n_30860_1_0 n_30870_1_0 n_30880_1_0 n_30890_1_0


foreach v of varlist n_30600_0_0-n_30890_0_0 {
	
format `v' %8.3f
}


/* MERGE POPULATION INCIDENCE RATES */

merge m:1 baseyear age sex using "CR_incid_mort.dta" 
drop if _merge==2 // ages outside the UKB age range
drop _merge


/* RUN THE CODE FOR THE DERIVED VARIABLES */

do "CR_read_data_derived.do"


/* EXCLUSIONS FOR ALL ANALYSES */

* exclude if outside 40-69 age range
drop if age<40 | age>=70 | age_calc<40 | age_calc>=70

* exclude if non-White British
*drop if ancestry_code!=9 // principal components ancestry

* exclude prevalent cases
drop if cr_timing==2

* exclude if no SNPs available
drop if snp_n==. 

* exclude people with Chron's disease
drop if ibd_chrons==1 

* exclude people with ulcerative colitis 
drop if ibd_colitis==1 

* exclude people with polyps
drop if polyp==1

* exclude if <6 weeks of follow-up
drop if agedeath!=. & deathdate<basedate+41 

* exclude if colorectal cancer in first six weeks
drop if cr==1 & cr_date<basedate+41


/* LIMIT TO UNRELATED INDIVIDUALS */

*** run this AFTER the exclusions so that we don't exclude people unecessarily

*** export the list of IDs to use in the R function ukb_gen_samples_to_remove
*** which is part of the ukbtools package
*** only run these lines if the list of people to exclude has been regenerated

*export delimited n_eid using "CR_related_check.csv", replace

*** then run the R script CR_related_R_IDs_to_drop.R 
*** this will generate the file CR_related_remove.csv

*import delimited "CR_related_remove.csv",  clear
*rename x n_eid
*save "CR_related_remove.dta", replace


merge 1:1 n_eid using "CR_related_remove.dta"
drop if _merge==3 | _merge==2
drop _merge


/* DIVIDE INTO TRAINING AND TESTING DATASETS */

splitsample if ancestry_code==9, generate(train_test) split(0.7 0.3) ///
	balance(cr sex) rseed(80171920)
recode train_test (.=2) // non white ancestries
label variable train_test "Study group"
label define train_test_label 1 "Training" 2 "Testing"  
label values train_test train_test_label


/* SAVE ANALYSIS FILE */

save "CR_analysis.dta", replace

