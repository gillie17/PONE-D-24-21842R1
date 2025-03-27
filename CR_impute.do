*** Colorectal cancer -- multiple imputation of missing data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

use "CR_analysis.dta", clear // analysis dataset

set maxiter 20 // don't want to wait for the default 300 iterations to fail

* training data
keep if train_test==1

* exclude if non-White British
drop if ancestry_code!=9 // principal components ancestry


drop incid_a0-mort_s_life84 // not needed in training data

* impute missing data
mi set mlong // use mi convert flong to use individual imputed datasets

mi register imputed ln_bmi ln_activity_mets walk_days walk_dur ///
	mod_days mod_dur vig_days vig_dur menop_hrt n_30690_0_0 n_30760_0_0 ///
	n_30780_0_0 n_30870_0_0 cr_deg1 screen10_new_time nsaid supp_vitamin_d ///
	supp_calcium alcohol smoke_ever fishoil diet_procmeat_cat ///
	diet_fruit_dried_cat diet_beef_cat diet_pork_cat diet_cereal ///
	diet_bread_whole diet_bread_white diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh diet_cereal_cat diet_bread_whole_cat ///
	diet_bread_white_cat diet_veg_cook_cat diet_veg_raw_cat diet_fruit_fresh_cat

mi register regular age_calc sex std_prs_xb diabetes screen10_new


*sample 10 // to do practice runs faster

display "$S_TIME  $S_DATE"

* 14 imputations 

mi impute chained (regress) ln_bmi ln_activity_mets n_30690_0_0 n_30760_0_0 ///
		n_30780_0_0 n_30870_0_0 ///
	(truncreg, ll(0) ul(10)) diet_fruit_fresh diet_veg_cook diet_veg_raw ///
	(truncreg, ll(0) ul(10) conditional(if screen10_new==1)) screen10_new_time /// 
	(pmm, knn(3)) diet_cereal diet_bread_whole diet_bread_white  ///
	(logit) cr_deg1 smoke_ever nsaid supp_vitamin_d supp_calcium fishoil ///
		diet_fruit_dried_cat ///
	(mlogit, conditional(if sex==0)) menop_hrt /// 
	(ologit) alcohol diet_procmeat_cat diet_beef_cat diet_pork_cat ///
	=cr age_calc std_prs_xb diabetes screen10_new, add(14) rseed(171920) ///
	dots augment force 
		
display "$S_TIME  $S_DATE"


* update the categorical diet variables

mi passive: generate diet_fruit_fresh_rd=round(diet_fruit_fresh) 
mi passive: generate diet_veg_cook_rd=round(diet_veg_cook) 
mi passive: generate diet_veg_raw_rd=round(diet_veg_raw) 

foreach n of numlist 1/14 {
	
replace diet_cereal_cat=0 if _mi_m==`n' & diet_cereal_cat==. & ///
	diet_cereal==0	
replace diet_cereal_cat=1 if _mi_m==`n' & diet_cereal_cat==. & ///
	inrange(diet_cereal,1,3)
replace diet_cereal_cat=2 if _mi_m==`n' & diet_cereal_cat==. & ///
	inrange(diet_cereal,4,6)
replace diet_cereal_cat=3 if _mi_m==`n' & diet_cereal_cat==. & ///
	inrange(diet_cereal,7,10)
		
replace diet_fruit_fresh_cat=0 if _mi_m==`n' & diet_fruit_fresh_cat==. & ///
	diet_fruit_fresh_rd==0
replace diet_fruit_fresh_cat=1 if _mi_m==`n' & diet_fruit_fresh_cat==. & ///
	diet_fruit_fresh_rd==1
replace diet_fruit_fresh_cat=2 if _mi_m==`n' & diet_fruit_fresh_cat==. & ///
	diet_fruit_fresh_rd==2
replace diet_fruit_fresh_cat=3 if _mi_m==`n' & diet_fruit_fresh_cat==. & ///
	inrange(diet_fruit_fresh_rd,3,10)		
		
replace diet_bread_whole_cat=0 if _mi_m==`n' & diet_bread_whole_cat==. & ///
	diet_bread_whole==0	
replace diet_bread_whole_cat=1 if _mi_m==`n' & diet_bread_whole_cat==. & ///
	inrange(diet_bread_whole,1,4)
replace diet_bread_whole_cat=2 if _mi_m==`n' & diet_bread_whole_cat==. & ///
	inrange(diet_bread_whole,5,10)
replace diet_bread_whole_cat=3 if _mi_m==`n' & diet_bread_whole_cat==. & ///
	inrange(diet_bread_whole,11,20)	
replace diet_bread_whole_cat=4 if _mi_m==`n' & diet_bread_whole_cat==. & ///
	inrange(diet_bread_whole,21,50)	
	
replace diet_bread_white_cat=0 if _mi_m==`n' & diet_bread_white_cat==. & ///
	diet_bread_white==0	
replace diet_bread_white_cat=1 if _mi_m==`n' & diet_bread_white_cat==. & ///
	inrange(diet_bread_white,1,4)
replace diet_bread_white_cat=2 if _mi_m==`n' & diet_bread_white_cat==. & ///
	inrange(diet_bread_white,5,10)
replace diet_bread_white_cat=3 if _mi_m==`n' & diet_bread_white_cat==. & ///
	inrange(diet_bread_white,11,20)	
replace diet_bread_white_cat=4 if _mi_m==`n' & diet_bread_white_cat==. & ///
	inrange(diet_bread_white,21,50)	

replace diet_veg_cook_cat=0 if _mi_m==`n' & diet_veg_cook_cat==. & ///
	diet_veg_cook_rd==0	
replace diet_veg_cook_cat=1 if _mi_m==`n' & diet_veg_cook_cat==. & ///
	diet_veg_cook_rd==1
replace diet_veg_cook_cat=2 if _mi_m==`n' & diet_veg_cook_cat==. & ///
	diet_veg_cook_rd==2
replace diet_veg_cook_cat=3 if _mi_m==`n' & diet_veg_cook_cat==. & ///
	inrange(diet_veg_cook_rd,3,10)

replace diet_veg_raw_cat=0 if _mi_m==`n' & diet_veg_raw_cat==. & ///
	diet_veg_raw_rd==0	
replace diet_veg_raw_cat=1 if _mi_m==`n' & diet_veg_raw_cat==. & ///
	diet_veg_raw_rd==1
replace diet_veg_raw_cat=2 if _mi_m==`n' & diet_veg_raw_cat==. & ///
	diet_veg_raw_rd==2
replace diet_veg_raw_cat=3 if _mi_m==`n' & diet_veg_raw_cat==. & ///
	inrange(diet_veg_raw_rd,3,10)
}


* new categories for bread

label define diet_bread_cat_new_label 0 "None" 1 "1-4 per week" ///
	2 "5-10 per week" 3 "11+ per week" 

mi passive: generate diet_bread_white_cat_new=diet_bread_white_cat
recode diet_bread_white_cat_new (4=3)
label variable diet_bread_white_cat_new "White bread slices per week: category"

mi passive: generate diet_bread_whole_cat_new=diet_bread_whole_cat
recode diet_bread_whole_cat_new (4=3)
label variable diet_bread_whole_cat_new "Wholemeal/grain bread slices per week: category"

label values diet_bread_white_cat_new diet_bread_whole_cat_new diet_bread_cat_new_label


* passive variables 

summ n_30760_0_0 if _mi_m==0
scalar sd30760=r(sd)
mi passive: generate n_30760_0_0_sd=n_30760_0_0/sd30760
label variable n_30760_0_0_sd "Per SD of HDL"

summ n_30870_0_0 if _mi_m==0
scalar sd30870=r(sd)
mi passive: generate n_30870_0_0_sd=n_30870_0_0/sd30870
label variable n_30870_0_0_sd "Per SD of tryglycerides"

summ n_30780_0_0 if _mi_m==0
scalar sd30780=r(sd)
mi passive: generate n_30780_0_0_sd=n_30780_0_0/sd30780
label variable n_30780_0_0_sd "Per SD of LDL"

summ n_30690_0_0 if _mi_m==0
scalar sd30690=r(sd)
mi passive: generate n_30690_0_0_sd=n_30690_0_0/sd30690
label variable n_30690_0_0_sd "Per SD of cholesterol"

summ ln_activity_mets if _mi_m==0
scalar sd_mets=r(sd)
mi passive: generate ln_activity_mets_sd=ln_activity_mets/sd_mets
label variable ln_activity_mets_sd "Per SD of activity METS"

summ ln_bmi if _mi_m==0
scalar sd_bmi=r(sd)
mi passive: generate ln_bmi_sd=ln_bmi/sd_bmi
label variable ln_bmi_sd "Per SD of BMI"


* centred values

summ ln_bmi if _mi_m==0
scalar ln_bmi_mean=r(mean)
mi passive: generate ln_bmi_c=ln_bmi-ln_bmi_mean
label variable ln_bmi_c "Centred ln_bmi"

summ age_calc if _mi_m==0
scalar age_calc_mean=r(mean)
mi passive: generate age_calc_c=age_calc-age_calc_mean
label variable age_calc_c "Centred age_calc"

summ n_30870_0_0 if _mi_m==0
scalar n_30870_0_0_mean=r(mean)
mi passive: generate n_30870_0_0_c=n_30870_0_0-n_30870_0_0_mean
label variable n_30870_0_0_c "Centred triglycerides"

summ n_30760_0_0 if _mi_m==0
scalar n_30760_0_0_mean=r(mean)
mi passive: generate n_30760_0_0_c=n_30760_0_0-n_30760_0_0_mean
label variable n_30760_0_0_c "Centred HDL"

summ n_30780_0_0 if _mi_m==0
scalar n_30780_0_0_mean=r(mean)
mi passive: generate n_30780_0_0_c=n_30780_0_0-n_30780_0_0_mean
label variable n_30780_0_0_c "Centred LDL"

summ n_30690_0_0 if _mi_m==0
scalar n_30690_0_0_mean=r(mean)
mi passive: generate n_30690_0_0_c=n_30690_0_0-n_30690_0_0_mean
label variable n_30690_0_0_c "Centred cholesterol"

summ ln_activity_mets if _mi_m==0
scalar ln_activity_mets_mean=r(mean)
mi passive: generate ln_activity_mets_c=ln_activity_mets-ln_activity_mets_mean
label variable ln_activity_mets_c "Centred activity"


* stset

mi stset cr_censor_age, id(n_eid) enter(time age_calc) failure(cr==1)


save "CR_analysis_imputed.dta", replace


