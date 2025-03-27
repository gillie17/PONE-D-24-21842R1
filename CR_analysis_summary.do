*** UK Biobank -- colorectal cancer summary analysis ***

version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Summary_`filedate'.log", name(summary) replace

use "CR_analysis.dta", clear // analysis dataset

drop if ancestry_code!=9 // UK ancestry only

tab train_test cr, miss co // all follow-up

* recodes

label define diet_bread_cat_new_label 0 "None" 1 "1-4 per week" ///
	2 "5-10 per week" 3 "11+ per week" 

generate diet_bread_white_cat_new=diet_bread_white_cat
recode diet_bread_white_cat_new (4=3)
label variable diet_bread_white_cat_new "White bread slices per week: category"

generate diet_bread_whole_cat_new=diet_bread_whole_cat
recode diet_bread_whole_cat_new (4=3)
label variable diet_bread_whole_cat_new "Wholemeal/grain bread slices per week: category"

label values diet_bread_white_cat_new diet_bread_whole_cat_new ///
	diet_bread_cat_new_label
	
	
*** number of SNPs genotyped

tab snp_n, miss

*** all follow-up -- right censored at 31/7/2019

tab cr sex, miss co ro

tabstat age_calc fup_time, column(statistics) labelwidth(20) varwidth(20) ///
	statistics(N mean sd min max) 

bysort cr: tabstat age_calc fup_time, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max) 

* age at diagnosis for all follow-up
tabstat cr_censor_age cr_fup_time if cr==1, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max) 

* female

tabstat age_calc fup_time if sex==0, column(statistics) labelwidth(20) ///
	varwidth(20) statistics(N mean sd min max) 

bysort cr: tabstat age_calc fup_time if sex==0, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max)

* age at diagnosis for all follow-up
tabstat cr_censor_age cr_fup_time if cr==1 & sex==0, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max)

* male

tabstat age_calc fup_time if sex==1, column(statistics) labelwidth(20) ///
	varwidth(20) statistics(N mean sd min max)

bysort cr: tabstat age_calc fup_time if sex==1, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max) 
	
* age at diagnosis for all follow-up
tabstat cr_censor_age cr_fup_time if cr==1 & sex==1, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max)
	
	
*** 10 years of follow-up

tab cr10 sex, miss co ro

tabstat age_calc fup10_time, column(statistics) varwidth(20) ///
	statistics(N mean sd min max) 

bysort cr10: tabstat age_calc fup10_time, column(statistics) varwidth(20) ///
	statistics(N mean sd min max) 

* age at diagnosis for 10 years of follow up	
tabstat cr_sir_censor10_age fup10_sir_time if cr10==1, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd min max) 

* female

tabstat age_calc fup10_time if sex==0, column(statistics) varwidth(20) ///
	statistics(N mean sd min max) 

bysort cr10: tabstat age_calc fup10_time if sex==0, column(statistics) ///
	varwidth(20) statistics(N mean sd min max) 

* age at diagnosis for 10 years of follow up	
tabstat cr_sir_censor10_age fup10_sir_time if cr10==1 & sex==0, ///
	column(statistics) labelwidth(20) varwidth(20) statistics(N mean sd min max) 	
	
* male

tabstat age_calc fup10_time if sex==1, column(statistics) varwidth(20) ///
	statistics(N mean sd min max) 

bysort cr10: tabstat age_calc fup10_time if sex==1, column(statistics) ///
	varwidth(20) statistics(N mean sd min max) 

* age at diagnosis for 10 years of follow up	
tabstat cr_sir_censor10_age fup10_sir_time if cr10==1 & sex==1, ///
	column(statistics) labelwidth(20) varwidth(20) statistics(N mean sd min max) 	

	
*** summary statistics

* all

* unaffected
misstable summarize sex lc_prs_xb ln_bmi activity_mets screen10_new_time ///
	n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh cr_deg1 screen10_new diabetes nsaid menop_hrt ///
	supp_calcium supp_vitamin_d fishoil alcohol smoke_ever diet_procmeat_cat ///
	diet_beef_cat diet_pork_cat diet_fruit_dried_cat diet_cereal_cat ///
	diet_bread_white_cat_new diet_bread_whole_cat_new if cr==0, all 

* affected
misstable summarize sex lc_prs_xb ln_bmi activity_mets screen10_new_time ///
	n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh cr_deg1 screen10_new diabetes nsaid menop_hrt ///
	supp_calcium supp_vitamin_d fishoil alcohol smoke_ever diet_procmeat_cat ///
	diet_beef_cat diet_pork_cat diet_fruit_dried_cat diet_cereal_cat ///
	diet_bread_white_cat_new diet_bread_whole_cat_new if cr==1, all 

tab sex cr, miss co
tab cr_deg1 cr, miss co

* missing values for mum, dad and sibling
bysort cr: tab miss_mum miss_dad if miss_sib==0, miss
bysort cr: tab miss_sib miss_mum if miss_dad==0, miss
bysort cr: tab miss_sib miss_dad if miss_mum==0, miss

tab screen10_new cr, miss co
tab diabetes cr, miss co
tab nsaid cr, miss co
tab menop_hrt cr if sex==0, miss co
tab supp_calcium cr, miss co
tab supp_vitamin_d cr, miss co
tab supp_fishoil cr, miss co
tab alcohol cr, miss co
tab smoke_ever cr, miss co
tab diet_procmeat_cat cr, miss co
tab diet_beef_cat cr, miss co
tab diet_pork_cat cr, miss co
tab diet_fruit_dried_cat cr, miss co
tab diet_cereal_cat cr, miss co
tab diet_bread_white_cat_new cr, miss co
tab diet_bread_whole_cat_new cr, miss co


bysort cr: tabstat screen10_new_time if screen10_new==1, column(statistics) ///
	labelwidth(20) varwidth(20) statistics(N mean sd p50 iqr min max) 

bysort cr: tabstat lc_prs_xb prs_rr bmi activity_mets n_30690_0_0 ///
	n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh, column(statistics) labelwidth(20) varwidth(20) ///
	statistics(N mean sd p50 iqr min max) 

	
* female

* unaffected
misstable summarize sex lc_prs_xb ln_bmi activity_mets screen10_new_time ///
	n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh cr_deg1 screen10_new diabetes nsaid menop_hrt ///
	supp_calcium supp_vitamin_d fishoil alcohol smoke_ever diet_procmeat_cat ///
	diet_beef_cat diet_pork_cat diet_fruit_dried_cat diet_cereal_cat ///
	diet_bread_white_cat_new diet_bread_whole_cat_new if cr==0 & sex==0, all 

* affected
misstable summarize sex lc_prs_xb ln_bmi activity_mets screen10_new_time ///
	n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh cr_deg1 screen10_new diabetes nsaid menop_hrt ///
	supp_calcium supp_vitamin_d fishoil alcohol smoke_ever diet_procmeat_cat ///
	diet_beef_cat diet_pork_cat diet_fruit_dried_cat diet_cereal_cat ///
	diet_bread_white_cat_new diet_bread_whole_cat_new if cr==1 & sex==0, all 

tab cr_deg1 cr if sex==0, miss co

* missing values for mum, dad and sibling
bysort cr: tab miss_mum miss_dad if miss_sib==0 & sex==0, miss
bysort cr: tab miss_sib miss_mum if miss_dad==0 & sex==0, miss
bysort cr: tab miss_sib miss_dad if miss_mum==0 & sex==0, miss

tab screen10_new cr if sex==0, miss co
tab diabetes cr if sex==0, miss co
tab nsaid cr if sex==0, miss co
tab menop_hrt cr if sex==0, miss co
tab supp_calcium cr if sex==0, miss co
tab supp_vitamin_d cr if sex==0, miss co
tab supp_fishoil cr if sex==0, miss co
tab alcohol cr if sex==0, miss co
tab smoke_ever cr if sex==0, miss co
tab diet_procmeat_cat cr if sex==0, miss co
tab diet_beef_cat cr if sex==0, miss co
tab diet_pork_cat cr if sex==0, miss co
tab diet_fruit_dried_cat cr if sex==0, miss co
tab diet_cereal_cat cr if sex==0, miss co
tab diet_bread_white_cat_new cr if sex==0, miss co
tab diet_bread_whole_cat_new cr if sex==0, miss co

bysort cr: tabstat screen10_new_time if sex==0 & screen10_new==1, ///
	column(statistics) labelwidth(20) varwidth(20) ///
	statistics(N mean sd p50 iqr min max) 

bysort cr: tabstat lc_prs_xb prs_rr bmi activity_mets n_30690_0_0 ///
	n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh if sex==0, column(statistics) labelwidth(20) ///
	varwidth(20) statistics(N mean sd p50 iqr min max) 

* male

* unaffected
misstable summarize sex lc_prs_xb ln_bmi activity_mets screen10_new_time ///
	n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh cr_deg1 screen10_new diabetes nsaid menop_hrt ///
	supp_calcium supp_vitamin_d fishoil alcohol smoke_ever diet_procmeat_cat ///
	diet_beef_cat diet_pork_cat diet_fruit_dried_cat diet_cereal_cat ///
	diet_bread_white_cat_new diet_bread_whole_cat_new if cr==0 & sex==1, all 

* affected
misstable summarize sex lc_prs_xb ln_bmi activity_mets screen10_new_time ///
	n_30690_0_0 n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh cr_deg1 screen10_new diabetes nsaid menop_hrt ///
	supp_calcium supp_vitamin_d fishoil alcohol smoke_ever diet_procmeat_cat ///
	diet_beef_cat diet_pork_cat diet_fruit_dried_cat diet_cereal_cat ///
	diet_bread_white_cat_new diet_bread_whole_cat_new if cr==1 & sex==1, all 

tab cr_deg1 cr if sex==1, miss co

* missing values for mum, dad and sibling
bysort cr: tab miss_mum miss_dad if miss_sib==0 & sex==1, miss
bysort cr: tab miss_sib miss_mum if miss_dad==0 & sex==1, miss
bysort cr: tab miss_sib miss_dad if miss_mum==0 & sex==1, miss

tab screen10_new cr if sex==1, miss co
tab diabetes cr if sex==1, miss co
tab nsaid cr if sex==1, miss co
tab supp_calcium cr if sex==1, miss co
tab supp_vitamin_d cr if sex==1, miss co
tab supp_fishoil cr if sex==1, miss co
tab alcohol cr if sex==1, miss co
tab smoke_ever cr if sex==1, miss co
tab diet_procmeat_cat cr if sex==1, miss co
tab diet_beef_cat cr if sex==1, miss co
tab diet_pork_cat cr if sex==1, miss co
tab diet_fruit_dried_cat cr if sex==1, miss co
tab diet_cereal_cat cr if sex==1, miss co
tab diet_bread_white_cat_new cr if sex==1, miss co
tab diet_bread_whole_cat_new cr if sex==1, miss co

bysort cr: tabstat screen10_new_time if sex==1 & screen10_new==1, ///
	column(statistics) labelwidth(20) varwidth(20) ///
	statistics(N mean sd p50 iqr min max) 

bysort cr: tabstat lc_prs_xb prs_rr bmi activity_mets n_30690_0_0 ///
	n_30760_0_0 n_30780_0_0 n_30870_0_0 diet_veg_cook diet_veg_raw ///
	diet_fruit_fresh if sex==1, column(statistics) labelwidth(20) ///
	varwidth(20) statistics(N mean sd p50 iqr min max) 

	
*** variables in models in testing dataset
	
misstable summarize std_prs_xb n_30870_0_0 if train_test==2 & sex==0, all 
misstable summarize std_prs_xb ln_bmi if train_test==2 & sex==1, all 
tab cr_deg1 sex if train_test==2, miss co
tab smoke_ever sex if train_test==2, miss co
tab screen10_new sex if train_test==2, miss co	
	
log close summary	
	
	

