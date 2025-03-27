version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

*log using "Results\Univariate.log", name(univariate) replace

use "CR_analysis.dta", clear // analysis dataset

* use all available follow-up for building the model in the training set

* training data
keep if train_test==1

* stset the data

stset cr_censor_age, id(n_eid) enter(time age_calc) failure(cr==1)


/* UNIVARIATE ANALYSES */

stcox sex // don't include in the model if no menop/hrt variable
estat phtest, detail
*estat phtest, plot(sex) name(sex, replace) // is okay

stcox std_prs_xb
estat phtest, detail

/*
estat phtest, plot(std_prs_xb) mcolor(eltblue) mlcolor(eltblue) msize(vsmall) ///
	lineopts(lcolor(ebblue) lwidth(medthick)) ///
	title("") xtitle("Age (years)", margin(t +2)) ///
	ytitle("Scaled Schoenfeld residuals:" "Standardised 140-SNP PRS") ///
	xlabel(40(5)80) yscale(range(-8 10)) ylabel(-8(2)10, angle(h)) ///
	note("") name(prs_xb, replace)
*/

stcox ln_bmi
*estat phtest, detail

*stcox exercise_m 
*estat phtest, detail

*stcox exercise_v
*estat phtest, detail

stcox i.cr_deg1_012
*estat phtest, detail

stcox activity_mets1000
*estat phtest, detail

stcox b2.walk_pace
*estat phtest, detail

stcox i.screen15_cat
*estat phtest, detail

stcox diabetes
*estat phtest, detail

stcox nsaid
*estat phtest, detail

stcox i.menop_hrt if sex==0
*estat phtest, detail

stcox supp_calcium
*estat phtest, detail

stcox supp_vitamin_d
*estat phtest, detail

stcox fishoil
*estat phtest, detail

stcox i.alcohol
*estat phtest, detail

stcox i.smoke_ever
*estat phtest, detail

stcox i.diet_procmeat_cat 
*estat phtest, detail

stcox i.diet_beef_cat 
*estat phtest, detail

stcox i.diet_pork_cat 
*estat phtest, detail

stcox i.diet_veg_cook_cat 
*estat phtest, detail

stcox i.diet_veg_raw_cat 
*estat phtest, detail

stcox i.diet_fruit_fresh_cat 
*estat phtest, detail

stcox diet_fruit_dried_cat 
*estat phtest, detail

stcox i.diet_cereal_cat 
*estat phtest, detail

stcox i.diet_bread_white_cat 
*estat phtest, detail

stcox i.diet_bread_whole_cat 
*estat phtest, detail

stcox n_30690_0_0 // cholesterol
*estat phtest, detail

stcox n_30760_0_0 // HDL
*estat phtest, detail

stcox n_30780_0_0 // LDL
*estat phtest, detail

stcox n_30870_0_0 // triglycerides
*estat phtest, detail

log close univariate
