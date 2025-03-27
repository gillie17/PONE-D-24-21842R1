version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Univariate_imputed_`filedate'.log", ///
	name(univariate_impute) replace

use "CR_analysis_imputed.dta", clear // imputed analysis dataset
mi convert flong // do this if doing separate analyses on the _mi_m


/* COMBINED */

mi estimate, hr: stcox std_prs_xb

mi estimate, hr: stcox ln_bmi_c

mi estimate, hr: stcox ln_activity_mets_c

mi estimate, hr: stcox n_30690_0_0_c // cholesterol

mi estimate, hr: stcox n_30760_0_0_c // HDL

mi estimate, hr: stcox n_30780_0_0_c // LDL

mi estimate, hr: stcox n_30870_0_0_c // triglycerides

mi estimate, hr: stcox diet_veg_cook

mi estimate, hr: stcox diet_veg_raw

mi estimate, hr: stcox diet_fruit_fresh

mi estimate, hr: stcox cr_deg1

mi estimate, hr: stcox screen10_new 

mi estimate, hr: stcox screen10_new_time if screen10_new==1

mi estimate, hr: stcox diabetes

mi estimate, hr: stcox nsaid

mi estimate, hr: stcox i.menop_hrt if sex==0

mi estimate, hr: stcox supp_calcium

mi estimate, hr: stcox supp_vitamin_d

mi estimate, hr: stcox fishoil

mi estimate, hr: stcox i.alcohol

mi estimate, hr: stcox smoke_ever

mi estimate, hr: stcox i.diet_procmeat_cat
mi test 1.diet_procmeat_cat 2.diet_procmeat_cat 3.diet_procmeat_cat

mi estimate, hr: stcox i.diet_beef_cat
mi test 1.diet_beef_cat 2.diet_beef_cat 3.diet_beef_cat

mi estimate, hr: stcox i.diet_pork_cat
mi test 1.diet_pork_cat 2.diet_pork_cat 3.diet_pork_cat

mi estimate, hr: stcox diet_fruit_dried_cat

mi estimate, hr: stcox i.diet_cereal_cat
mi test 1.diet_cereal_cat 2.diet_cereal_cat 3.diet_cereal_cat

mi estimate, hr: stcox i.diet_bread_white_cat_new
mi test 1.diet_bread_white_cat 2.diet_bread_white_cat 3.diet_bread_white_cat 
	
mi estimate, hr: stcox i.diet_bread_whole_cat_new
mi test 1.diet_bread_whole_cat 2.diet_bread_whole_cat 3.diet_bread_whole_cat 


/* WOMEN */

mi estimate, hr: stcox std_prs_xb if sex==0

mi estimate, hr: stcox ln_bmi_c if sex==0

mi estimate, hr: stcox ln_activity_mets_c if sex==0

mi estimate, hr: stcox n_30690_0_0_c if sex==0 // cholesterol

mi estimate, hr: stcox n_30760_0_0_c if sex==0 // HDL

mi estimate, hr: stcox n_30780_0_0_c if sex==0 // LDL

mi estimate, hr: stcox n_30870_0_0_c if sex==0 // triglycerides

mi estimate, hr: stcox diet_veg_cook if sex==0

mi estimate, hr: stcox diet_veg_raw if sex==0

mi estimate, hr: stcox diet_fruit_fresh if sex==0

mi estimate, hr: stcox cr_deg1 if sex==0

mi estimate, hr: stcox screen10_new if sex==0

mi estimate, hr: stcox screen10_new_time if sex==0 & screen10_new==1

mi estimate, hr: stcox diabetes if sex==0

mi estimate, hr: stcox nsaid if sex==0

mi estimate, hr: stcox i.menop_hrt if sex==0

mi estimate, hr: stcox supp_calcium if sex==0

mi estimate, hr: stcox supp_vitamin_d if sex==0

mi estimate, hr: stcox fishoil if sex==0

mi estimate, hr: stcox i.alcohol if sex==0

mi estimate, hr: stcox smoke_ever if sex==0

mi estimate, hr: stcox i.diet_procmeat_cat if sex==0
mi test 1.diet_procmeat_cat 2.diet_procmeat_cat 3.diet_procmeat_cat

mi estimate, hr: stcox i.diet_beef_cat if sex==0
mi test 1.diet_beef_cat 2.diet_beef_cat 3.diet_beef_cat

mi estimate, hr: stcox i.diet_pork_cat if sex==0
mi test 1.diet_pork_cat 2.diet_pork_cat 3.diet_pork_cat

mi estimate, hr: stcox diet_fruit_dried_cat if sex==0

mi estimate, hr: stcox i.diet_cereal_cat if sex==0
mi test 1.diet_cereal_cat 2.diet_cereal_cat 3.diet_cereal_cat

mi estimate, hr: stcox i.diet_bread_white_cat_new if sex==0
mi test 1.diet_bread_white_cat 2.diet_bread_white_cat 3.diet_bread_white_cat 
	
mi estimate, hr: stcox i.diet_bread_whole_cat_new if sex==0
mi test 1.diet_bread_whole_cat 2.diet_bread_whole_cat 3.diet_bread_whole_cat 


/* MEN */

mi estimate, hr: stcox std_prs_xb if sex==1

mi estimate, hr: stcox ln_bmi_c if sex==1

mi estimate, hr: stcox ln_activity_mets_c if sex==1

mi estimate, hr: stcox n_30690_0_0_c if sex==1 // cholesterol

mi estimate, hr: stcox n_30760_0_0_c if sex==1 // HDL

mi estimate, hr: stcox n_30780_0_0_c if sex==1 // LDL

mi estimate, hr: stcox n_30870_0_0_c if sex==1 // triglycerides

mi estimate, hr: stcox diet_veg_cook if sex==1

mi estimate, hr: stcox diet_veg_raw if sex==1

mi estimate, hr: stcox diet_fruit_fresh if sex==1

mi estimate, hr: stcox cr_deg1 if sex==1

mi estimate, hr: stcox screen10_new if sex==1

mi estimate, hr: stcox screen10_new_time if sex==1 & screen10_new==1

mi estimate, hr: stcox diabetes if sex==1

mi estimate, hr: stcox nsaid if sex==1

mi estimate, hr: stcox supp_calcium if sex==1

mi estimate, hr: stcox supp_vitamin_d if sex==1

mi estimate, hr: stcox fishoil if sex==1

mi estimate, hr: stcox i.alcohol if sex==1

mi estimate, hr: stcox smoke_ever if sex==1

mi estimate, hr: stcox i.diet_procmeat_cat if sex==1
mi test 1.diet_procmeat_cat 2.diet_procmeat_cat 3.diet_procmeat_cat

mi estimate, hr: stcox i.diet_beef_cat if sex==1
mi test 1.diet_beef_cat 2.diet_beef_cat 3.diet_beef_cat

mi estimate, hr: stcox i.diet_pork_cat if sex==1
mi test 1.diet_pork_cat 2.diet_pork_cat 3.diet_pork_cat

mi estimate, hr: stcox diet_fruit_dried_cat if sex==1

mi estimate, hr: stcox i.diet_cereal_cat if sex==1
mi test 1.diet_cereal_cat 2.diet_cereal_cat 3.diet_cereal_cat

mi estimate, hr: stcox i.diet_bread_white_cat_new if sex==1
mi test 1.diet_bread_white_cat 2.diet_bread_white_cat 3.diet_bread_white_cat 
	
mi estimate, hr: stcox i.diet_bread_whole_cat_new if sex==1
mi test 1.diet_bread_whole_cat 2.diet_bread_whole_cat 3.diet_bread_whole_cat 

log close univariate_impute