* how many imputations to do

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\HowManyImputations_`filedate'.log", ///
	name(how_many) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m


* number of imputations needed - calculated after the first 11 imputations

* all variables, no interaction
mi estimate, imputations(1/11): stcox std_prs_xb sex cr_deg1 ln_bmi diabetes ///
	i.alcohol smoke_ever ln_activity_mets i.menop_hrt supp_vitamin_d ///
	supp_calcium nsaid screen10_new screen10_new_time i.diet_procmeat_cat ///
	i.diet_beef_cat i.diet_pork_cat i.diet_cereal_cat i.diet_bread_whole_cat ///
	i.diet_bread_white_cat diet_veg_cook diet_veg_raw diet_fruit_fresh ///
	diet_fruit_dried_cat n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0, ///
	baselevels hr 

how_many_imputations, cv_se(0.05) 

* omitting activity
mi estimate, imputations(1/11): stcox std_prs_xb sex cr_deg1 ln_bmi diabetes i.alcohol ///
	smoke_ever i.menop_hrt supp_vitamin_d supp_calcium nsaid ///
	screen10_new screen10_new_time i.diet_procmeat_cat i.diet_beef_cat ///
	i.diet_pork_cat i.diet_cereal_cat i.diet_bread_whole_cat ///
	i.diet_bread_white_cat diet_veg_cook diet_veg_raw diet_fruit_fresh ///
	diet_fruit_dried_cat n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0, ///
	baselevels hr 

how_many_imputations, cv_se(0.05) 

* omitting blood lipids
mi estimate, imputations(1/11): stcox std_prs_xb sex cr_deg1 ln_bmi diabetes i.alcohol ///
	smoke_ever i.menop_hrt supp_vitamin_d supp_calcium nsaid ///
	screen10_new screen10_new_time i.diet_procmeat_cat i.diet_beef_cat ///
	i.diet_pork_cat i.diet_cereal_cat i.diet_bread_whole_cat ///
	i.diet_bread_white_cat diet_veg_cook diet_veg_raw diet_fruit_fresh ///
	diet_fruit_dried_cat, baselevels hr 

how_many_imputations, cv_se(0.05) 

* women - final model

mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	screen10_new i.menop_hrt if sex==0

how_many_imputations, cv_se(0.05) 

* men - final model

mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat screen10_new ///
	screen10_new_time if sex==1

how_many_imputations, cv_se(0.05) 


log close how_many
