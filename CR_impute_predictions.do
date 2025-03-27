
*** UK Biobank -- colorectal cancer predictions ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Final_predictions_`filedate'.log", ///
	name(final_predict) replace

use "CR_analysis_imputed.dta", clear
mi convert flong // do this if doing predictions

* predicted values from model estimates
mi predict final_f_xb using models\model_kevin_f if sex==0, xb
mi predict final_m_xb using models\model_kevin_m if sex==1, xb


*** female - training data ***

estimates use models\model_kevin_f
stcox, nohr

generate final_f_calc = _b[std_prs_xb]*std_prs_xb + ///
	_b[cr_deg1]*cr_deg1 + ///
	_b[smoke_ever]*smoke_ever + ///
	_b[screen10_new]*screen10_new + ///
	_b[n_30870_0_0_c]*n_30870_0_0_c if sex==0 & _mi_m==0

capture program drop female_betas
	
program female_betas

di as text _dup(40) "-"
di "lnHR coefficients"
di as text _dup(40) "-"
di "                std_prs_xb = "_b[std_prs_xb]
di "                   cr_deg1 = "_b[cr_deg1]
di "                smoke_ever = "_b[smoke_ever]
di "              screen15_cat = "_b[screen10_new]
di "             n_30870_0_0_c = "_b[n_30870_0_0_c]

di as text _dup(40) "-"

end

female_betas	
	
* predictions are in _mi_m==0		
mi xeq 0: summarize final_f_xb final_f_calc if sex==0 & final_f_calc!=.
	
	
*** male - training data ***

estimates use models\model_kevin_m
stcox, nohr

generate final_m_calc = _b[std_prs_xb]*std_prs_xb + ///
	_b[cr_deg1]*cr_deg1 + ///
	_b[smoke_ever]*smoke_ever + ///
	_b[screen10_new]*screen10_new + ///
	_b[ln_bmi_c]*ln_bmi_c if sex==1 & _mi_m==0


capture program drop male_betas	
	
program male_betas

di as text _dup(40) "-"
di "lnHR coefficients"
di as text _dup(40) "-"
di "                std_prs_xb = "_b[std_prs_xb]
di "                   cr_deg1 = "_b[cr_deg1]
di "                smoke_ever = "_b[smoke_ever]
di "              screen15_cat = "_b[screen10_new]
di "                  ln_bmi_c = "_b[ln_bmi_c]
di as text _dup(40) "-"

end

male_betas
	
* predictions are in _mi_m==0		
mi xeq 0: summarize final_m_xb final_m_calc if sex==1 & final_m_calc!=.


log close final_predict

