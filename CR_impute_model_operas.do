*** UK Biobank -- colorectal cancer operas using imputed data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Model_opera_imputed_`filedate'.log", ///
	name(model_opera_impute) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m

mi query // counts the imputations
local M = r(M)

foreach N of numlist 1/`M' {

*** women ***

* FH model
regress std_prs_xb cr_deg1 if sex==0 & _mi_m==`N'
predict fh_prs_resid_f`N', residuals 

logit cr_deg1 std_prs_xb if sex==0 & _mi_m==`N'
predict fh_deg1_resid_f`N', residuals 

summ fh_prs_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_fh_prs_resid_f`N'=r(sd)
generate fh_prs_opera_f`N'=fh_prs_resid_f`N'/sd_fh_prs_resid_f`N' ///
	if sex==0 & _mi_m==`N'

summ fh_deg1_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_fh_deg1_resid_f`N'=r(sd)
generate fh_deg1_opera_f`N'=fh_deg1_resid_f`N'/sd_fh_deg1_resid_f`N' ///
	if sex==0 & _mi_m==`N'

* multivariable model
regress std_prs_xb cr_deg1 smoke_ever screen10_new n_30870_0_0_c ///
	if sex==0 & _mi_m==`N'
predict prs_resid_f`N', residuals 

logit cr_deg1 std_prs_xb smoke_ever screen10_new n_30870_0_0_c ///
	if sex==0 & _mi_m==`N'
predict deg1_resid_f`N', residuals 

logit smoke_ever std_prs_xb cr_deg1 screen10_new n_30870_0_0_c ///
	if sex==0 & _mi_m==`N'
predict smoke_resid_f`N', residuals 

logit screen10_new std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	if sex==0 & _mi_m==`N'
predict screen_resid_f`N', residuals 

regress n_30870_0_0_c std_prs_xb cr_deg1 smoke_ever screen10_new ///
	if sex==0 & _mi_m==`N'
predict trigly_resid_f`N', residuals 

summ prs_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_prs_resid_f`N'=r(sd)
generate prs_opera_f`N'=prs_resid_f`N'/sd_prs_resid_f`N' ///
	if sex==0 & _mi_m==`N'

summ deg1_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_deg1_resid_f`N'=r(sd)
generate deg1_opera_f`N'=deg1_resid_f`N'/sd_deg1_resid_f`N' ///
	if sex==0 & _mi_m==`N'

summ smoke_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_smoke_resid_f`N'=r(sd)
generate smoke_opera_f`N'=smoke_resid_f`N'/sd_smoke_resid_f`N' ///
	if sex==0 & _mi_m==`N'

summ screen_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_screen_resid_f`N'=r(sd)
generate screen_opera_f`N'=screen_resid_f`N'/sd_screen_resid_f`N' ///
	if sex==0 & _mi_m==`N'

summ trigly_resid_f`N' if sex==0 & _mi_m==`N'
scalar sd_trigly_resid_f`N'=r(sd)
generate trigly_opera_f`N'=trigly_resid_f`N'/sd_trigly_resid_f`N' ///
	if sex==0 & _mi_m==`N'


*** men ***

* FH model
regress std_prs_xb cr_deg1 if sex==1 & _mi_m==`N'
predict fh_prs_resid_m`N', residuals 

logit cr_deg1 std_prs_xb if sex==1 & _mi_m==`N'
predict fh_deg1_resid_m`N', residuals 

summ fh_prs_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_mh_prs_resid_m`N'=r(sd)
generate fh_prs_opera_m`N'=fh_prs_resid_m`N'/sd_mh_prs_resid_m`N' ///
	if sex==1 & _mi_m==`N'

summ fh_deg1_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_mh_deg1_resid_m`N'=r(sd)
generate fh_deg1_opera_m`N'=fh_deg1_resid_m`N'/sd_mh_deg1_resid_m`N' ///
	if sex==1 & _mi_m==`N'

* multivariable model
regress std_prs_xb cr_deg1 smoke_ever screen10_new ln_bmi_c ///
	if sex==1 & _mi_m==`N'
predict prs_resid_m`N', residuals 

logit cr_deg1 std_prs_xb smoke_ever screen10_new ln_bmi_c ///
	if sex==1 & _mi_m==`N'
predict deg1_resid_m`N', residuals 

logit smoke_ever std_prs_xb cr_deg1 screen10_new ln_bmi_c ///
	if sex==1 & _mi_m==`N'
predict smoke_resid_m`N', residuals 

logit screen10_new std_prs_xb cr_deg1 smoke_ever ln_bmi_c ///
	if sex==1 & _mi_m==`N'
predict screen_resid_m`N', residuals 

regress ln_bmi_c std_prs_xb cr_deg1 smoke_ever screen10_new ///
	if sex==1 & _mi_m==`N'
predict bmi_resid_m`N', residuals 

summ prs_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_prs_resid_m`N'=r(sd)
generate prs_opera_m`N'=prs_resid_m`N'/sd_prs_resid_m`N' ///
	if sex==1 & _mi_m==`N'

summ deg1_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_deg1_resid_m`N'=r(sd)
generate deg1_opera_m`N'=deg1_resid_m`N'/sd_deg1_resid_m`N' ///
	if sex==1 & _mi_m==`N'

summ smoke_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_smoke_resid_m`N'=r(sd)
generate smoke_opera_m`N'=smoke_resid_m`N'/sd_smoke_resid_m`N' ///
	if sex==1 & _mi_m==`N'

summ screen_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_screen_resid_m`N'=r(sd)
generate screen_opera_m`N'=screen_resid_m`N'/sd_screen_resid_m`N' ///
	if sex==1 & _mi_m==`N'

summ bmi_resid_m`N' if sex==1 & _mi_m==`N'
scalar sd_bmi_resid_m`N'=r(sd)
generate bmi_opera_m`N'=bmi_resid_m`N'/sd_bmi_resid_m`N' ///
	if sex==1 & _mi_m==`N'	
}	
	
	
**# OPERAS ***

* FAMILY HISTORY MODEL - women

scalar fh_prs_beta_f=0 // start with 0
scalar fh_deg1_beta_f=0

scalar fh_prs_sesq_f=0	
scalar fh_deg1_sesq_f=0		
	
foreach N of numlist 1/`M' {

* women
stcox fh_prs_opera_f`N' fh_deg1_opera_f`N' if sex==0 // model

scalar fh_prs_beta_f`N'=_b[fh_prs_opera_f`N'] // estimate
scalar fh_deg1_beta_f`N'=_b[fh_deg1_opera_f`N'] 

scalar fh_prs_se_f`N'=_se[fh_prs_opera_f`N'] // standard error
scalar fh_deg1_se_f`N'=_se[fh_deg1_opera_f`N'] 

scalar fh_prs_beta_f = fh_prs_beta_f + fh_prs_beta_f`N' // estimate sum
scalar fh_deg1_beta_f = fh_deg1_beta_f + fh_deg1_beta_f`N' 

scalar fh_prs_sesq_f = fh_prs_sesq_f + (fh_prs_se_f`N')^2 // variance sum
scalar fh_deg1_sesq_f = fh_deg1_sesq_f + (fh_deg1_se_f`N')^2 

}

scalar fh_prs_beta_f_pooled = fh_prs_beta_f / `M' // pooled estimate
scalar fh_deg1_beta_f_pooled = fh_deg1_beta_f / `M' 

scalar fh_prs_beta_f_diff=0 // (estimate - pooled estimate)^2
scalar fh_deg1_beta_f_diff=0 


foreach N of numlist 1/`M' {
	
scalar fh_prs_beta_f_diff = fh_prs_beta_f_diff + ///
	(fh_prs_beta_f`N' - fh_prs_beta_f_pooled)^2
	
scalar fh_deg1_beta_f_diff = fh_deg1_beta_f_diff + ///
	(fh_deg1_beta_f`N' - fh_deg1_beta_f_pooled)^2

}

scalar vw_fh_prs_f = fh_prs_sesq_f/`M'
scalar vb_fh_prs_f = fh_prs_beta_f_diff/(`M'-1)

scalar vw_fh_deg1_f = fh_deg1_sesq_f/`M'
scalar vb_fh_deg1_f = fh_deg1_beta_f_diff/(`M'-1)

scalar fh_prs_se_f_pooled = sqrt(vw_fh_prs_f + vb_fh_prs_f + vb_fh_prs_f/`M')
scalar fh_deg1_se_f_pooled = sqrt(vw_fh_deg1_f + vb_fh_deg1_f + vb_fh_deg1_f/`M')


* FULL MODEL - women

scalar prs_beta_f=0 // start with 0
scalar deg1_beta_f=0
scalar smoke_beta_f=0
scalar screen_beta_f=0
scalar trigly_beta_f=0

scalar prs_sesq_f=0	
scalar deg1_sesq_f=0
scalar smoke_sesq_f=0
scalar screen_sesq_f=0
scalar trigly_sesq_f=0	
	
	
foreach N of numlist 1/`M' {

* women
stcox prs_opera_f`N' deg1_opera_f`N' smoke_opera_f`N' screen_opera_f`N' ///
	trigly_opera_f`N' if sex==0 // model

scalar prs_beta_f`N'=_b[prs_opera_f`N'] // estimate
scalar deg1_beta_f`N'=_b[deg1_opera_f`N'] 
scalar smoke_beta_f`N'=_b[smoke_opera_f`N'] 
scalar screen_beta_f`N'=_b[screen_opera_f`N'] 
scalar trigly_beta_f`N'=_b[trigly_opera_f`N'] 

scalar prs_se_f`N'=_se[prs_opera_f`N'] // standard error
scalar deg1_se_f`N'=_se[deg1_opera_f`N'] 
scalar smoke_se_f`N'=_se[smoke_opera_f`N'] 
scalar screen_se_f`N'=_se[screen_opera_f`N'] 
scalar trigly_se_f`N'=_se[trigly_opera_f`N'] 

scalar prs_beta_f = prs_beta_f + prs_beta_f`N' // estimate sum
scalar deg1_beta_f = deg1_beta_f + deg1_beta_f`N' 
scalar smoke_beta_f = smoke_beta_f + smoke_beta_f`N' 
scalar screen_beta_f = screen_beta_f + screen_beta_f`N' 
scalar trigly_beta_f = trigly_beta_f + trigly_beta_f`N' 

scalar prs_sesq_f = prs_sesq_f + (prs_se_f`N')^2 // variance sum
scalar deg1_sesq_f = deg1_sesq_f + (deg1_se_f`N')^2  
scalar smoke_sesq_f = smoke_sesq_f + (smoke_se_f`N')^2  
scalar screen_sesq_f = screen_sesq_f + (screen_se_f`N')^2  
scalar trigly_sesq_f = trigly_sesq_f + (trigly_se_f`N')^2  

}

scalar prs_beta_f_pooled = prs_beta_f / `M' // pooled estimate
scalar deg1_beta_f_pooled = deg1_beta_f / `M' 
scalar smoke_beta_f_pooled = smoke_beta_f / `M' 
scalar screen_beta_f_pooled = screen_beta_f / `M' 
scalar trigly_beta_f_pooled = trigly_beta_f / `M' 

scalar prs_beta_f_diff=0 // (estimate - pooled estimate)^2
scalar deg1_beta_f_diff=0 
scalar smoke_beta_f_diff=0 
scalar screen_beta_f_diff=0 
scalar trigly_beta_f_diff=0 

foreach N of numlist 1/`M' {
	
scalar prs_beta_f_diff = prs_beta_f_diff + ///
	(prs_beta_f`N' - prs_beta_f_pooled)^2
	
scalar deg1_beta_f_diff = deg1_beta_f_diff + ///
	(deg1_beta_f`N' - deg1_beta_f_pooled)^2

scalar smoke_beta_f_diff = smoke_beta_f_diff + ///
	(smoke_beta_f`N' - smoke_beta_f_pooled)^2
	
scalar screen_beta_f_diff = screen_beta_f_diff + ///
	(screen_beta_f`N' - screen_beta_f_pooled)^2
	
scalar trigly_beta_f_diff = trigly_beta_f_diff + ///
	(trigly_beta_f`N' - trigly_beta_f_pooled)^2
	
}

scalar vw_prs_f = prs_sesq_f/`M'
scalar vw_deg1_f = deg1_sesq_f/`M'
scalar vw_smoke_f = smoke_sesq_f/`M'
scalar vw_screen_f = screen_sesq_f/`M'
scalar vw_trigly_f = trigly_sesq_f/`M'

scalar vb_prs_f = prs_beta_f_diff/(`M'-1)
scalar vb_deg1_f = deg1_beta_f_diff/(`M'-1)
scalar vb_smoke_f = smoke_beta_f_diff/(`M'-1)
scalar vb_screen_f = screen_beta_f_diff/(`M'-1)
scalar vb_trigly_f = trigly_beta_f_diff/(`M'-1)


scalar prs_se_f_pooled = sqrt(vw_prs_f + vb_prs_f + vb_prs_f/`M')
scalar deg1_se_f_pooled = sqrt(vw_deg1_f + vb_deg1_f + vb_deg1_f/`M')
scalar smoke_se_f_pooled = sqrt(vw_smoke_f + vb_smoke_f + vb_smoke_f/`M')
scalar screen_se_f_pooled = sqrt(vw_screen_f + vb_screen_f + vb_screen_f/`M')
scalar trigly_se_f_pooled = sqrt(vw_trigly_f + vb_trigly_f + vb_trigly_f/`M')

**# FAMILY HISTORY MODEL - men

scalar fh_prs_beta_m=0 // start with 0
scalar fh_deg1_beta_m=0

scalar fh_prs_sesq_m=0	
scalar fh_deg1_sesq_m=0		
	
foreach N of numlist 1/`M' {

* men
stcox fh_prs_opera_m`N' fh_deg1_opera_m`N' if sex==1 // model

scalar fh_prs_beta_m`N'=_b[fh_prs_opera_m`N'] // estimate
scalar fh_deg1_beta_m`N'=_b[fh_deg1_opera_m`N'] 

scalar fh_prs_se_m`N'=_se[fh_prs_opera_m`N'] // standard error
scalar fh_deg1_se_m`N'=_se[fh_deg1_opera_m`N'] 

scalar fh_prs_beta_m = fh_prs_beta_m + fh_prs_beta_m`N' // estimate sum
scalar fh_deg1_beta_m = fh_deg1_beta_m + fh_deg1_beta_m`N' 

scalar fh_prs_sesq_m = fh_prs_sesq_m + (fh_prs_se_m`N')^2 // variance sum
scalar fh_deg1_sesq_m = fh_deg1_sesq_m + (fh_deg1_se_m`N')^2 

}

scalar fh_prs_beta_m_pooled = fh_prs_beta_m / `M' // pooled estimate
scalar fh_deg1_beta_m_pooled = fh_deg1_beta_m / `M' 

scalar fh_prs_beta_m_diff=0 // (estimate - pooled estimate)^2
scalar fh_deg1_beta_m_diff=0 


foreach N of numlist 1/`M' {
	
scalar fh_prs_beta_m_diff = fh_prs_beta_m_diff + ///
	(fh_prs_beta_m`N' - fh_prs_beta_m_pooled)^2
	
scalar fh_deg1_beta_m_diff = fh_deg1_beta_m_diff + ///
	(fh_deg1_beta_m`N' - fh_deg1_beta_m_pooled)^2

}

scalar vw_fh_prs_m = fh_prs_sesq_m/`M'
scalar vb_fh_prs_m = fh_prs_beta_m_diff/(`M'-1)

scalar vw_fh_deg1_m = fh_deg1_sesq_m/`M'
scalar vb_fh_deg1_m = fh_deg1_beta_m_diff/(`M'-1)

scalar fh_prs_se_m_pooled = sqrt(vw_fh_prs_m + vb_fh_prs_m + vb_fh_prs_m/`M')
scalar fh_deg1_se_m_pooled = sqrt(vw_fh_deg1_m + vb_fh_deg1_m + vb_fh_deg1_m/`M')


**# FULL MODEL - men

scalar prs_beta_m=0 // start with 0
scalar deg1_beta_m=0
scalar smoke_beta_m=0
scalar screen_beta_m=0
scalar bmi_beta_m=0

scalar prs_sesq_m=0	
scalar deg1_sesq_m=0
scalar smoke_sesq_m=0
scalar screen_sesq_m=0
scalar bmi_sesq_m=0	
	
	
foreach N of numlist 1/`M' {

* men
stcox prs_opera_m`N' deg1_opera_m`N' smoke_opera_m`N' screen_opera_m`N' ///
	bmi_opera_m`N' if sex==1 // model

scalar prs_beta_m`N'=_b[prs_opera_m`N'] // estimate
scalar deg1_beta_m`N'=_b[deg1_opera_m`N'] 
scalar smoke_beta_m`N'=_b[smoke_opera_m`N'] 
scalar screen_beta_m`N'=_b[screen_opera_m`N'] 
scalar bmi_beta_m`N'=_b[bmi_opera_m`N'] 

scalar prs_se_m`N'=_se[prs_opera_m`N'] // standard error
scalar deg1_se_m`N'=_se[deg1_opera_m`N'] 
scalar smoke_se_m`N'=_se[smoke_opera_m`N'] 
scalar screen_se_m`N'=_se[screen_opera_m`N'] 
scalar bmi_se_m`N'=_se[bmi_opera_m`N'] 

scalar prs_beta_m = prs_beta_m + prs_beta_m`N' // estimate sum
scalar deg1_beta_m = deg1_beta_m + deg1_beta_m`N' 
scalar smoke_beta_m = smoke_beta_m + smoke_beta_m`N' 
scalar screen_beta_m = screen_beta_m + screen_beta_m`N' 
scalar bmi_beta_m = bmi_beta_m + bmi_beta_m`N' 

scalar prs_sesq_m = prs_sesq_m + (prs_se_m`N')^2 // variance sum
scalar deg1_sesq_m = deg1_sesq_m + (deg1_se_m`N')^2  
scalar smoke_sesq_m = smoke_sesq_m + (smoke_se_m`N')^2  
scalar screen_sesq_m = screen_sesq_m + (screen_se_m`N')^2  
scalar bmi_sesq_m = bmi_sesq_m + (bmi_se_m`N')^2  

}

scalar prs_beta_m_pooled = prs_beta_m / `M' // pooled estimate
scalar deg1_beta_m_pooled = deg1_beta_m / `M' 
scalar smoke_beta_m_pooled = smoke_beta_m / `M' 
scalar screen_beta_m_pooled = screen_beta_m / `M' 
scalar bmi_beta_m_pooled = bmi_beta_m / `M' 

scalar prs_beta_m_diff=0 // (estimate - pooled estimate)^2
scalar deg1_beta_m_diff=0 
scalar smoke_beta_m_diff=0 
scalar screen_beta_m_diff=0 
scalar bmi_beta_m_diff=0 

foreach N of numlist 1/`M' {
	
scalar prs_beta_m_diff = prs_beta_m_diff + ///
	(prs_beta_m`N' - prs_beta_m_pooled)^2
	
scalar deg1_beta_m_diff = deg1_beta_m_diff + ///
	(deg1_beta_m`N' - deg1_beta_m_pooled)^2

scalar smoke_beta_m_diff = smoke_beta_m_diff + ///
	(smoke_beta_m`N' - smoke_beta_m_pooled)^2
	
scalar screen_beta_m_diff = screen_beta_m_diff + ///
	(screen_beta_m`N' - screen_beta_m_pooled)^2
	
scalar bmi_beta_m_diff = bmi_beta_m_diff + ///
	(bmi_beta_m`N' - bmi_beta_m_pooled)^2
	
}

scalar vw_prs_m = prs_sesq_m/`M'
scalar vw_deg1_m = deg1_sesq_m/`M'
scalar vw_smoke_m = smoke_sesq_m/`M'
scalar vw_screen_m = screen_sesq_m/`M'
scalar vw_bmi_m = bmi_sesq_m/`M'

scalar vb_prs_m = prs_beta_m_diff/(`M'-1)
scalar vb_deg1_m = deg1_beta_m_diff/(`M'-1)
scalar vb_smoke_m = smoke_beta_m_diff/(`M'-1)
scalar vb_screen_m = screen_beta_m_diff/(`M'-1)
scalar vb_bmi_m = bmi_beta_m_diff/(`M'-1)


scalar prs_se_m_pooled = sqrt(vw_prs_m + vb_prs_m + vb_prs_m/`M')
scalar deg1_se_m_pooled = sqrt(vw_deg1_m + vb_deg1_m + vb_deg1_m/`M')
scalar smoke_se_m_pooled = sqrt(vw_smoke_m + vb_smoke_m + vb_smoke_m/`M')
scalar screen_se_m_pooled = sqrt(vw_screen_m + vb_screen_m + vb_screen_m/`M')
scalar bmi_se_m_pooled = sqrt(vw_bmi_m + vb_bmi_m + vb_bmi_m/`M')


**# DISPLAY RESULTS

capture program drop trick

program trick
di 
di "Family history model - women"
di "----------------------------"
di
di "FH prs_opera_f - beta      " round(fh_prs_beta_f_pooled,0.001)
di "FH prs_opera_f - se        " round(fh_prs_se_f_pooled,0.001) 

di 
di "FH prs_opera_f - HR        " round(exp(fh_prs_beta_f_pooled),0.001)
di "FH prs_opera_f - 95% CI    " round(exp(fh_prs_beta_f_pooled - ///
	1.96*fh_prs_se_f_pooled),0.001) " to " round(exp(fh_prs_beta_f_pooled + ///
	1.96*fh_prs_se_f_pooled),0.001)
di "FH prs_opera_f - P-value   " round(2*(1-normal(fh_prs_beta_f_pooled / ///
	fh_prs_se_f_pooled)),0.001)
di
di "FH deg1_opera_f - beta:    " round(fh_deg1_beta_f_pooled,0.001)
di "FH deg1_opera_f = se:      " round(fh_deg1_se_f_pooled,0.001)
di
di "FH deg1_opera_f - HR        " round(exp(fh_deg1_beta_f_pooled),0.001)
di "FH deg1_opera_f - 95% CI    " round(exp(fh_deg1_beta_f_pooled - ///
	1.96*fh_deg1_se_f_pooled),0.001) " to " round(exp(fh_deg1_beta_f_pooled + ///
	1.96*fh_deg1_se_f_pooled),0.001)
di "FH deg1_opera_f - P-value   " round(2*(1-normal(fh_deg1_beta_f_pooled / ///
	fh_deg1_se_f_pooled)),0.001)
di 
di "Full model - women"
di "------------------"
di
di "FULL prs_opera_f - beta      " round(prs_beta_f_pooled,0.001)
di "FULL prs_opera_f - se        " round(prs_se_f_pooled,0.001) 

di 
di "FULL prs_opera_f - HR        " round(exp(prs_beta_f_pooled),0.001)
di "FULL prs_opera_f - 95% CI    " round(exp(prs_beta_f_pooled - ///
	1.96*prs_se_f_pooled),0.001) " to " round(exp(prs_beta_f_pooled + ///
	1.96*prs_se_f_pooled),0.001)
di "FULL prs_opera_f - P-value   " round(2*(1-normal(prs_beta_f_pooled / ///
	prs_se_f_pooled)),0.001)
di
di "FULL deg1_opera_f - beta:    " round(deg1_beta_f_pooled,0.001)
di "FULL deg1_opera_f = se:      " round(deg1_se_f_pooled,0.001)
di
di "FULL deg1_opera_f - HR        " round(exp(deg1_beta_f_pooled),0.001)
di "FULL deg1_opera_f - 95% CI    " round(exp(deg1_beta_f_pooled - ///
	1.96*deg1_se_f_pooled),0.001) " to " round(exp(deg1_beta_f_pooled + ///
	1.96*deg1_se_f_pooled),0.001)
di "FULL deg1_opera_f - P-value   " round(2*(1-normal(deg1_beta_f_pooled / ///
	deg1_se_f_pooled)),0.001)
di
di "FULL smoke_opera_f - beta:    " round(smoke_beta_f_pooled,0.001)
di "FULL smoke_opera_f = se:      " round(smoke_se_f_pooled,0.001)
di
di "FULL smoke_opera_f - HR        " round(exp(smoke_beta_f_pooled),0.001)
di "FULL smoke_opera_f - 95% CI    " round(exp(smoke_beta_f_pooled - ///
	1.96*smoke_se_f_pooled),0.001) " to " round(exp(smoke_beta_f_pooled + ///
	1.96*smoke_se_f_pooled),0.001)
di "FULL smoke_opera_f - P-value   " round(2*(1-normal(smoke_beta_f_pooled / ///
	smoke_se_f_pooled)),0.001)
di
di "FULL screen_opera_f - beta:    " round(screen_beta_f_pooled,0.001)
di "FULL screen_opera_f = se:      " round(screen_se_f_pooled,0.001)
di
di "FULL screen_opera_f - HR        " round(exp(screen_beta_f_pooled),0.001)
di "FULL screen_opera_f - 95% CI    " round(exp(screen_beta_f_pooled - ///
	1.96*screen_se_f_pooled),0.001) " to " round(exp(screen_beta_f_pooled + ///
	1.96*screen_se_f_pooled),0.001)
di "FULL screen_opera_f - P-value   " round(2*(normal(screen_beta_f_pooled / ///
	screen_se_f_pooled)),0.001)
di
di "FULL trigly_opera_f - beta:    " round(trigly_beta_f_pooled,0.001)
di "FULL trigly_opera_f = se:      " round(trigly_se_f_pooled,0.001)
di
di "FULL trigly_opera_f - HR        " round(exp(trigly_beta_f_pooled),0.001)
di "FULL trigly_opera_f - 95% CI    " round(exp(trigly_beta_f_pooled - ///
	1.96*trigly_se_f_pooled),0.001) " to " round(exp(trigly_beta_f_pooled + ///
	1.96*trigly_se_f_pooled),0.001)
di "FULL trigly_opera_f - P-value   " round(2*(1-normal(trigly_beta_f_pooled / ///
	trigly_se_f_pooled)),0.001)
di 
di "Family history model - men"
di "----------------------------"
di
di "FH prs_opera_m - beta      " round(fh_prs_beta_m_pooled,0.001)
di "FH prs_opera_m - se        " round(fh_prs_se_m_pooled,0.001) 

di 
di "FH prs_opera_m - HR        " round(exp(fh_prs_beta_m_pooled),0.001)
di "FH prs_opera_m - 95% CI    " round(exp(fh_prs_beta_m_pooled - ///
	1.96*fh_prs_se_m_pooled),0.001) " to " round(exp(fh_prs_beta_m_pooled + ///
	1.96*fh_prs_se_m_pooled),0.001)
di "FH prs_opera_m - P-value   " round(2*(1-normal(fh_prs_beta_m_pooled / ///
	fh_prs_se_m_pooled)),0.001)
di
di "FH deg1_opera_m - beta:    " round(fh_deg1_beta_m_pooled,0.001)
di "FH deg1_opera_m = se:      " round(fh_deg1_se_m_pooled,0.001)
di
di "FH deg1_opera_m - HR        " round(exp(fh_deg1_beta_m_pooled),0.001)
di "FH deg1_opera_m - 95% CI    " round(exp(fh_deg1_beta_m_pooled - ///
	1.96*fh_deg1_se_m_pooled),0.001) " to " round(exp(fh_deg1_beta_m_pooled + ///
	1.96*fh_deg1_se_m_pooled),0.001)
di "FH deg1_opera_m - P-value   " round(2*(1-normal(fh_deg1_beta_m_pooled / ///
	fh_deg1_se_m_pooled)),0.001)
di 
di "Full model - men"
di "------------------"
di
di "FULL prs_opera_m - beta      " round(prs_beta_m_pooled,0.001)
di "FULL prs_opera_m - se        " round(prs_se_m_pooled,0.001) 

di 
di "FULL prs_opera_m - HR        " round(exp(prs_beta_m_pooled),0.001)
di "FULL prs_opera_m - 95% CI    " round(exp(prs_beta_m_pooled - ///
	1.96*prs_se_m_pooled),0.001) " to " round(exp(prs_beta_m_pooled + ///
	1.96*prs_se_m_pooled),0.001)
di "FULL prs_opera_m - P-value   " round(2*(1-normal(prs_beta_m_pooled / ///
	prs_se_m_pooled)),0.001)
di
di "FULL deg1_opera_m - beta:    " round(deg1_beta_m_pooled,0.001)
di "FULL deg1_opera_m = se:      " round(deg1_se_m_pooled,0.001)
di
di "FULL deg1_opera_m - HR        " round(exp(deg1_beta_m_pooled),0.001)
di "FULL deg1_opera_m - 95% CI    " round(exp(deg1_beta_m_pooled - ///
	1.96*deg1_se_m_pooled),0.001) " to " round(exp(deg1_beta_m_pooled + ///
	1.96*deg1_se_m_pooled),0.001)
di "FULL deg1_opera_m - P-value   " round(2*(1-normal(deg1_beta_m_pooled / ///
	deg1_se_m_pooled)),0.001)
di
di "FULL smoke_opera_m - beta:    " round(smoke_beta_m_pooled,0.001)
di "FULL smoke_opera_m = se:      " round(smoke_se_m_pooled,0.001)
di
di "FULL smoke_opera_m - HR        " round(exp(smoke_beta_m_pooled),0.001)
di "FULL smoke_opera_m - 95% CI    " round(exp(smoke_beta_m_pooled - ///
	1.96*smoke_se_m_pooled),0.001) " to " round(exp(smoke_beta_m_pooled + ///
	1.96*smoke_se_m_pooled),0.001)
di "FULL smoke_opera_m - P-value   " round(2*(1-normal(smoke_beta_m_pooled / ///
	smoke_se_m_pooled)),0.001)
di
di "FULL screen_opera_m - beta:    " round(screen_beta_m_pooled,0.001)
di "FULL screen_opera_m = se:      " round(screen_se_m_pooled,0.001)
di
di "FULL screen_opera_m - HR        " round(exp(screen_beta_m_pooled),0.001)
di "FULL screen_opera_m - 95% CI    " round(exp(screen_beta_m_pooled - ///
	1.96*screen_se_m_pooled),0.001) " to " round(exp(screen_beta_m_pooled + ///
	1.96*screen_se_m_pooled),0.001)
di "FULL screen_opera_m - P-value   " round(2*(normal(screen_beta_m_pooled / ///
	screen_se_m_pooled)),0.001)
di
di "FULL bmi_opera_m - beta:    " round(bmi_beta_m_pooled,0.001)
di "FULL bmi_opera_m = se:      " round(bmi_se_m_pooled,0.001)
di
di "FULL bmi_opera_m - HR        " round(exp(bmi_beta_m_pooled),0.001)
di "FULL bmi_opera_m - 95% CI    " round(exp(bmi_beta_m_pooled - ///
	1.96*bmi_se_m_pooled),0.001) " to " round(exp(bmi_beta_m_pooled + ///
	1.96*bmi_se_m_pooled),0.001)
di "FULL bmi_opera_m - P-value   " round(2*(1-normal(bmi_beta_m_pooled / ///
	bmi_se_m_pooled)),0.001)	
	
end

trick

log close model_opera_impute

exit



