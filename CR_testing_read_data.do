*** colorectal cancer -- testing analysis ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

use "CR_analysis.dta", clear // analysis dataset

*** centred values ***

* BMI
summ ln_bmi
scalar ln_bmi_mean=r(mean)
generate ln_bmi_c=ln_bmi-ln_bmi_mean
label variable ln_bmi_c "Centred ln_bmi"

* triglycerides
summ n_30870_0_0
scalar n_30870_0_0_mean=r(mean)
generate n_30870_0_0_c=n_30870_0_0-n_30870_0_0_mean
label variable n_30870_0_0_c "Centred triglycerides"

* age at baseline
summ age_calc if train_test==1 // use the same value as in the training set
scalar age_calc_mean=r(mean)  
generate age_calc_c=age_calc-age_calc_mean
label variable age_calc_c "Centred age_calc"


*** missing data ***

misstable summarize std_prs_xb cr_deg1 smoke_ever screen10_new n_30870_0_0_c ///
	if sex==0, all
	
misstable summarize std_prs_xb cr_deg1 smoke_ever screen10_new ln_bmi_c ///
	if sex==1, all

*** deal with missing data ***

recode cr_deg1 (.=0) 
recode smoke_ever (.=0) 
recode n_30870_0_0_c (.=0) 
recode ln_bmi_c (.=0) 


*** ancestry groups ***

generate ancestry_group=.
replace ancestry_group=1 if ancestry_code==9 // UK
replace ancestry_group=2 if inlist(ancestry_code,1,5,6,8) // other Europe
replace ancestry_group=3 if inlist(ancestry_code,2,7) // Black
replace ancestry_group=4 if ancestry_code==4 // South Asian
replace ancestry_group=5 if ancestry_code==3 // East Asian
label variable ancestry_group "Collapsed ancestry groups"
label define ancestry_group_label 1 "UK" 2 "other Europe" 3 "Black" ///
	4 "South Asian" 5 "East Asian"
label values ancestry_group ancestry_group_label

*** limit to testing dataset ***

drop if train_test==1
	
	
*** Aviv's family history and PRS model - RR calculation ***

generate crc_fh_prs_rr=prs_rr * 0.92 if cr_deg1==0 | cr_deg1==.
replace crc_fh_prs_rr=prs_rr * 2.10 if cr_deg1==1


*** Aviv's family history - RR calculation ***

generate crc_fh_rr=0.92 if cr_deg1==0 | cr_deg1==.
replace crc_fh_rr=2.10 if cr_deg1==1


*** Aviv's PRS - RR calculation ***

generate crc_prs_rr=prs_rr 


*** new FH and PRS model - RR calculation ***

* women
estimates use models\model_fh_prs_f // use the saved estimates
stcox // to replay the estimates
predict newfhprs_f_xb if sex==0, xb

summ newfhprs_f_xb if sex==0 & ancestry_code==9
scalar newfhprs_f_xb_mean = r(mean)

generate crc_newfh_prs_rr = exp(newfhprs_f_xb - newfhprs_f_xb_mean) if sex==0 
	
* men
estimates use models\model_fh_prs_m // use the saved estimates
stcox // to replay the estimates
predict newfhprs_m_xb if sex==1, xb

summ newfhprs_m_xb if sex==1 & ancestry_code==9
scalar newfhprs_m_xb_mean = r(mean)

replace crc_newfh_prs_rr = exp(newfhprs_m_xb - newfhprs_m_xb_mean) if sex==1 


*** new FH model - RR calculation ***

* women
estimates use models\model_fh_f // use the saved estimates
stcox // to replay the estimates
predict newfh_f_xb if sex==0, xb

summ newfh_f_xb if sex==0 & ancestry_code==9
scalar newfh_f_xb_mean = r(mean)

generate crc_newfh_rr = exp(newfh_f_xb - newfh_f_xb_mean) if sex==0 
	
* men
estimates use models\model_fh_m // use the saved estimates
stcox // to replay the estimates
predict newfh_m_xb if sex==1, xb

summ newfh_m_xb if sex==1 & ancestry_code==9
scalar newfh_m_xb_mean = r(mean)

replace crc_newfh_rr = exp(newfh_m_xb - newfh_m_xb_mean) if sex==1 


*** new PRS model - RR calculation ***

* women
estimates use models\model_prs_f // use the saved estimates
stcox // to replay the estimates
predict newprs_f_xb if sex==0, xb

summ newprs_f_xb if sex==0 & ancestry_code==9
scalar newprs_f_xb_mean = r(mean)

generate crc_newprs_rr = exp(newprs_f_xb - newprs_f_xb_mean) if sex==0 
	
* men
estimates use models\model_prs_m // use the saved estimates
stcox // to replay the estimates
predict newprs_m_xb if sex==1, xb

summ newprs_m_xb if sex==1 & ancestry_code==9
scalar newprs_m_xb_mean = r(mean)

replace crc_newprs_rr = exp(newprs_m_xb - newprs_m_xb_mean) if sex==1 


*** new model - RR calculation ***

* women
estimates use models\model_kevin_f // use the saved estimates
stcox // to replay the estimates
predict final_f_xb if sex==0, xb

summ final_f_xb if sex==0 & ancestry_code==9
scalar final_f_xb_mean = r(mean)

generate crc_new_rr = exp(final_f_xb - final_f_xb_mean) if sex==0 

* men
estimates use models\model_kevin_m // use the saved estimates
stcox // to replay the estimates
predict final_m_xb if sex==1, xb

summ final_m_xb if sex==1 & ancestry_code==9
scalar final_m_xb_mean = r(mean)

replace crc_new_rr = exp(final_m_xb - final_m_xb_mean) if sex==1 



/* ABSOLUTE 10-YEAR RISK SCORES -- with competing mortality adjustment */


*** average risk (no RR in equation) ***

foreach n of numlist 0/9 {

generate crc_av_10_a`n' = ((incid_a`n') / (incid_a`n' + mort_a`n')) * ///
	(exp(-incid_s_a`n') / exp(-incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-incid_a`n' - mort_a`n'))
}

egen crc_av_10yr=rowtotal(crc_av_10_a0 - crc_av_10_a9) 


*** Aviv's family history and PRS model ***

foreach n of numlist 0/9 {

generate crc_fh_prs_10_a`n' = ((crc_fh_prs_rr*incid_a`n') / ///
		(crc_fh_prs_rr*incid_a`n' + mort_a`n')) * ///
	(exp(-crc_fh_prs_rr*incid_s_a`n') / exp(-crc_fh_prs_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_fh_prs_rr*incid_a`n' - mort_a`n'))
}

egen crc_fh_prs_10yr=rowtotal(crc_fh_prs_10_a0 - crc_fh_prs_10_a9) 


*** Aviv's family history model ***

foreach n of numlist 0/9 {

generate crc_fh_10_a`n' = ((crc_fh_rr*incid_a`n') / ///
		(crc_fh_rr*incid_a`n' + mort_a`n')) * ///
	(exp(-crc_fh_rr*incid_s_a`n') / exp(-crc_fh_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_fh_rr*incid_a`n' - mort_a`n'))
}

egen crc_fh_10yr=rowtotal(crc_fh_10_a0 - crc_fh_10_a9) 


*** Aviv's PRS model ***

foreach n of numlist 0/9 {

generate crc_prs_10_a`n' = ((crc_prs_rr*incid_a`n') / ///
		(crc_prs_rr*incid_a`n' + mort_a`n')) * ///
	(exp(-crc_prs_rr*incid_s_a`n') / exp(-crc_prs_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_prs_rr*incid_a`n' - mort_a`n'))
}

egen crc_prs_10yr=rowtotal(crc_prs_10_a0 - crc_prs_10_a9) 

*** new family history and PRS model ***

foreach n of numlist 0/9 {

generate crc_newfh_prs_10_a`n' = ((crc_newfh_prs_rr*incid_a`n') / ///
		(crc_newfh_prs_rr*incid_a`n' + mort_a`n')) * ///
	(exp(-crc_newfh_prs_rr*incid_s_a`n') / ///
	exp(-crc_newfh_prs_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_newfh_prs_rr*incid_a`n' - mort_a`n'))
}

egen crc_newfh_prs_10yr=rowtotal(crc_newfh_prs_10_a0 - crc_newfh_prs_10_a9) 


*** new family history model ***

foreach n of numlist 0/9 {

generate crc_newfh_10_a`n' = ((crc_newfh_rr*incid_a`n') / ///
		(crc_newfh_rr*incid_a`n' + mort_a`n')) * ///
	(exp(-crc_newfh_rr*incid_s_a`n') / ///
	exp(-crc_newfh_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_newfh_rr*incid_a`n' - mort_a`n'))
}

egen crc_newfh_10yr=rowtotal(crc_newfh_10_a0 - crc_newfh_10_a9) 


*** new PRS model ***

foreach n of numlist 0/9 {

generate crc_newprs_10_a`n' = ((crc_newprs_rr*incid_a`n') / ///
		(crc_newprs_rr*incid_a`n' + mort_a`n')) * ///
	(exp(-crc_newprs_rr*incid_s_a`n') / ///
	exp(-crc_newprs_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_newprs_rr*incid_a`n' - mort_a`n'))
}

egen crc_newprs_10yr=rowtotal(crc_newprs_10_a0 - crc_newprs_10_a9) 


*** new model ***

foreach n of numlist 0/9 {

generate crc_new_10_a`n' = ((crc_new_rr*incid_a`n') / (crc_new_rr*incid_a`n' ///
		+ mort_a`n')) * ///
	(exp(-crc_new_rr*incid_s_a`n') / exp(-crc_new_rr*incid_s_a0)) * ///
	(exp(-mort_s_a`n') / exp(-mort_s_a0)) * ///
	(1-exp(-crc_new_rr*incid_a`n' - mort_a`n'))
}

egen crc_new_10yr=rowtotal(crc_new_10_a0 - crc_new_10_a9) 



/* FULL-LIFETIME RISK SCORES -- with competing mortality adjustment */

*** average risk (no RR in equation)

foreach n of numlist 0/84 {

generate crc_av_life`n' = ((incid_life`n') / ///
		(incid_life`n' + mort_life`n')) * ///
	(exp(-incid_s_life`n') / exp(-incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-incid_life`n' - mort_life`n'))
}

egen crc_av_full=rowtotal(crc_av_life0 - crc_av_life84) 


*** Aviv's family history and PRS model ***
 
foreach n of numlist 0/84 {

generate crc_fh_prs_life`n' = ((crc_fh_prs_rr*incid_life`n') / ///
		(crc_fh_prs_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_fh_prs_rr*incid_s_life`n') / ///
	exp(-crc_fh_prs_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_fh_prs_rr*incid_life`n' - mort_life`n'))
}

egen crc_fh_prs_full=rowtotal(crc_fh_prs_life0 - crc_fh_prs_life84) 


*** Aviv's family history model ***
 
foreach n of numlist 0/84 {

generate crc_fh_life`n' = ((crc_fh_rr*incid_life`n') / ///
		(crc_fh_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_fh_rr*incid_s_life`n') / ///
	exp(-crc_fh_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_fh_rr*incid_life`n' - mort_life`n'))
}

egen crc_fh_full=rowtotal(crc_fh_life0 - crc_fh_life84) 


*** Aviv's PRS model ***
 
foreach n of numlist 0/84 {

generate crc_prs_life`n' = ((crc_prs_rr*incid_life`n') / ///
		(crc_prs_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_prs_rr*incid_s_life`n') / ///
	exp(-crc_prs_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_prs_rr*incid_life`n' - mort_life`n'))
}

egen crc_prs_full=rowtotal(crc_prs_life0 - crc_prs_life84) 



*** new family history and PRS model ***
 
foreach n of numlist 0/84 {

generate crc_newfh_prs_life`n' = ((crc_newfh_prs_rr*incid_life`n') / ///
		(crc_newfh_prs_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_newfh_prs_rr*incid_s_life`n') / ///
	exp(-crc_newfh_prs_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_newfh_prs_rr*incid_life`n' - mort_life`n'))
}

egen crc_newfh_prs_full=rowtotal(crc_newfh_prs_life0 - crc_newfh_prs_life84) 


*** new family history model ***
 
foreach n of numlist 0/84 {

generate crc_newfh_life`n' = ((crc_newfh_rr*incid_life`n') / ///
		(crc_newfh_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_newfh_rr*incid_s_life`n') / ///
	exp(-crc_newfh_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_newfh_rr*incid_life`n' - mort_life`n'))
}

egen crc_newfh_full=rowtotal(crc_newfh_life0 - crc_newfh_life84) 


*** new PRS model ***
 
foreach n of numlist 0/84 {

generate crc_newprs_life`n' = ((crc_newprs_rr*incid_life`n') / ///
		(crc_newprs_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_newprs_rr*incid_s_life`n') / ///
	exp(-crc_newprs_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_newprs_rr*incid_life`n' - mort_life`n'))
}

egen crc_newprs_full=rowtotal(crc_newprs_life0 - crc_newprs_life84) 


*** new model ***
 
foreach n of numlist 0/84 {

generate crc_new_life`n' = ((crc_new_rr*incid_life`n') / ///
		(crc_new_rr*incid_life`n' + mort_life`n')) * ///
	(exp(-crc_new_rr*incid_s_life`n') / exp(-crc_new_rr*incid_s_life0)) * ///
	(exp(-mort_s_life`n') / exp(-mort_s_life0)) * ///
	(1-exp(-crc_new_rr*incid_life`n' - mort_life`n'))
}

egen crc_new_full=rowtotal(crc_new_life0 - crc_new_life84) 



/* TURN ABSOLUTE RISKS INTO LOG ODDS */

* 10-year risk
generate crc_av_10yr_logit=logit(crc_av_10yr)
generate crc_fh_10yr_logit=logit(crc_fh_10yr)
generate crc_prs_10yr_logit=logit(crc_prs_10yr)
generate crc_fh_prs_10yr_logit=logit(crc_fh_prs_10yr)
generate crc_newfh_10yr_logit=logit(crc_newfh_10yr)
generate crc_newprs_10yr_logit=logit(crc_newprs_10yr)
generate crc_newfh_prs_10yr_logit=logit(crc_newfh_prs_10yr)
generate crc_new_10yr_logit=logit(crc_new_10yr)


generate crc_av_full_logit=logit(crc_av_full)
generate crc_fh_full_logit=logit(crc_fh_full)
generate crc_prs_full_logit=logit(crc_prs_full)
generate crc_fh_prs_full_logit=logit(crc_fh_prs_full)
generate crc_newfh_full_logit=logit(crc_newfh_full)
generate crc_newprs_full_logit=logit(crc_newprs_full)
generate crc_newfh_prs_full_logit=logit(crc_newfh_prs_full)
generate crc_new_full_logit=logit(crc_new_full)


/* RISK PER SD */

* Average risk
summarize crc_av_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_av_10yr_logit=r(sd)
generate crc_av_10yr_logit_sd=crc_av_10yr_logit/sdcrc_av_10yr_logit if sex==0

summarize crc_av_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_av_10yr_logit=r(sd)
replace crc_av_10yr_logit_sd=crc_av_10yr_logit/sdcrc_av_10yr_logit if sex==1

* Aviv's FH and PRS
summarize crc_fh_prs_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_fh_prs_10yr_logit=r(sd)
generate crc_fh_prs_10yr_logit_sd=crc_fh_prs_10yr_logit/sdcrc_fh_prs_10yr_logit if sex==0

summarize crc_fh_prs_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_fh_prs_10yr_logit=r(sd)
replace crc_fh_prs_10yr_logit_sd=crc_fh_prs_10yr_logit/sdcrc_fh_prs_10yr_logit if sex==1

* Aviv's FH
summarize crc_fh_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_fh_10yr_logit=r(sd)
generate crc_fh_10yr_logit_sd=crc_fh_10yr_logit/sdcrc_fh_10yr_logit if sex==0

summarize crc_fh_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_fh_10yr_logit=r(sd)
replace crc_fh_10yr_logit_sd=crc_fh_10yr_logit/sdcrc_fh_10yr_logit if sex==1

* Aviv's PRS
summarize crc_prs_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_prs_10yr_logit=r(sd)
generate crc_prs_10yr_logit_sd=crc_prs_10yr_logit/sdcrc_prs_10yr_logit if sex==0

summarize crc_prs_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_prs_10yr_logit=r(sd)
replace crc_prs_10yr_logit_sd=crc_prs_10yr_logit/sdcrc_prs_10yr_logit if sex==1

* new FH and PRS
summarize crc_newfh_prs_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_newfh_prs_10yr_logit=r(sd)
generate crc_newfh_prs_10yr_logit_sd=crc_newfh_prs_10yr_logit/sdcrc_newfh_prs_10yr_logit if sex==0

summarize crc_newfh_prs_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_newfh_prs_10yr_logit=r(sd)
replace crc_newfh_prs_10yr_logit_sd=crc_newfh_prs_10yr_logit/sdcrc_newfh_prs_10yr_logit if sex==1

* new FH 
summarize crc_newfh_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_newfh_10yr_logit=r(sd)
generate crc_newfh_10yr_logit_sd=crc_newfh_10yr_logit/sdcrc_newfh_10yr_logit if sex==0

summarize crc_newfh_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_newfh_10yr_logit=r(sd)
replace crc_newfh_10yr_logit_sd=crc_newfh_10yr_logit/sdcrc_newfh_10yr_logit if sex==1

* new PRS
summarize crc_newprs_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_newprs_10yr_logit=r(sd)
generate crc_newprs_10yr_logit_sd=crc_newprs_10yr_logit/sdcrc_newprs_10yr_logit if sex==0

summarize crc_newprs_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_newprs_10yr_logit=r(sd)
replace crc_newprs_10yr_logit_sd=crc_newprs_10yr_logit/sdcrc_newprs_10yr_logit if sex==1

* new model
summarize crc_new_10yr_logit if sex==0 & ancestry_code==9
scalar sdcrc_new_10yr_logit=r(sd)
generate crc_new_10yr_logit_sd=crc_new_10yr_logit/sdcrc_new_10yr_logit if sex==0

summarize crc_new_10yr_logit if sex==1 & ancestry_code==9
scalar sdcrc_new_10yr_logit=r(sd)
replace crc_new_10yr_logit_sd=crc_new_10yr_logit/sdcrc_new_10yr_logit if sex==1


*** quintiles of risk ***

* women

pctile crc_av_10yr_quintile_f = crc_av_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_av_10yr_quint_f = crc_av_10yr, ///
	cutpoints(crc_av_10yr_quintile_f)
generate crc_av_10yr_quint=crc_av_10yr_quint_f if sex==0
label variable crc_av_10yr_quint "Average 10-year quintiles"

pctile crc_fh_prs_10yr_quintile_f = crc_fh_prs_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_fh_prs_10yr_quint_f = crc_fh_prs_10yr, ///
	cutpoints(crc_fh_prs_10yr_quintile_f)
generate crc_fh_prs_10yr_quint=crc_fh_prs_10yr_quint_f if sex==0
label variable crc_fh_prs_10yr_quint "FH & PRS model 10-year quintiles"

pctile crc_fh_10yr_quintile_f = crc_fh_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_fh_10yr_quint_f = crc_fh_10yr, ///
	cutpoints(crc_fh_10yr_quintile_f)
generate crc_fh_10yr_quint=crc_fh_10yr_quint_f if sex==0
label variable crc_fh_10yr_quint "FH model 10-year quintiles"

pctile crc_prs_10yr_quintile_f = crc_prs_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_prs_10yr_quint_f = crc_prs_10yr, ///
	cutpoints(crc_prs_10yr_quintile_f)
generate crc_prs_10yr_quint=crc_prs_10yr_quint_f if sex==0
label variable crc_prs_10yr_quint "PRS model 10-year quintiles"

pctile crc_newfh_prs_10yr_quintile_f = crc_newfh_prs_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_newfh_prs_10yr_quint_f = crc_newfh_prs_10yr, ///
	cutpoints(crc_newfh_prs_10yr_quintile_f)
generate crc_newfh_prs_10yr_quint=crc_newfh_prs_10yr_quint_f if sex==0
label variable crc_newfh_prs_10yr_quint "New FH & PRS model 10-year quintiles"

pctile crc_newfh_10yr_quintile_f = crc_newfh_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_newfh_10yr_quint_f = crc_newfh_10yr, cutpoints(crc_newfh_10yr_quintile_f)
generate crc_newfh_10yr_quint=crc_newfh_10yr_quint_f if sex==0
label variable crc_newfh_10yr_quint "New FH model 10-year quintiles"

pctile crc_newprs_10yr_quintile_f = crc_newprs_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_newprs_10yr_quint_f = crc_newprs_10yr, ///
	cutpoints(crc_newprs_10yr_quintile_f)
generate crc_newprs_10yr_quint=crc_newprs_10yr_quint_f if sex==0
label variable crc_newfh_prs_10yr_quint "New PRS model 10-year quintiles"

pctile crc_new_10yr_quintile_f = crc_new_10yr if sex==0 & ///
	ancestry_code==9, nq(5)
xtile crc_new_10yr_quint_f = crc_new_10yr, cutpoints(crc_new_10yr_quintile_f)
generate crc_new_10yr_quint=crc_new_10yr_quint_f if sex==0
label variable crc_new_10yr_quint "New model 10-year quintiles"

* men
pctile crc_av_10yr_quintile_m = crc_av_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_av_10yr_quint_m = crc_av_10yr, ///
	cutpoints(crc_av_10yr_quintile_m)
replace crc_av_10yr_quint=crc_av_10yr_quint_m if sex==1

pctile crc_fh_prs_10yr_quintile_m = crc_fh_prs_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_fh_prs_10yr_quint_m = crc_fh_prs_10yr, ///
	cutpoints(crc_fh_prs_10yr_quintile_m)
replace crc_fh_prs_10yr_quint=crc_fh_prs_10yr_quint_m if sex==1

pctile crc_fh_10yr_quintile_m = crc_fh_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_fh_10yr_quint_m = crc_fh_10yr, ///
	cutpoints(crc_fh_10yr_quintile_m)
replace crc_fh_10yr_quint=crc_fh_10yr_quint_m if sex==1

pctile crc_prs_10yr_quintile_m = crc_prs_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_prs_10yr_quint_m = crc_prs_10yr, ///
	cutpoints(crc_prs_10yr_quintile_m)
replace crc_prs_10yr_quint=crc_prs_10yr_quint_m if sex==1

pctile crc_newfh_prs_10yr_quintile_m = crc_newfh_prs_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_newfh_prs_10yr_quint_m = crc_newfh_prs_10yr, ///
	cutpoints(crc_newfh_prs_10yr_quintile_m)
replace crc_newfh_prs_10yr_quint=crc_newfh_prs_10yr_quint_m if sex==1

pctile crc_newfh_10yr_quintile_m = crc_newfh_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_newfh_10yr_quint_m = crc_newfh_10yr, cutpoints(crc_newfh_10yr_quintile_m)
replace crc_newfh_10yr_quint=crc_newfh_10yr_quint_m if sex==1

pctile crc_newprs_10yr_quintile_m = crc_newprs_10yr if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_newprs_10yr_quint_m = crc_newprs_10yr, ///
	cutpoints(crc_newprs_10yr_quintile_m)
replace crc_newprs_10yr_quint=crc_newprs_10yr_quint_m if sex==1

pctile crc_new_10yr_quintile_m = crc_new_10yr  if sex==1 & ///
	ancestry_code==9, nq(5)
xtile crc_new_10yr_quint_m = crc_new_10yr, cutpoints(crc_new_10yr_quintile_m)
replace crc_new_10yr_quint=crc_new_10yr_quint_m if sex==1

drop crc_av_10yr_quintile_f crc_av_10yr_quint_f ///
	crc_fh_prs_10yr_quintile_f crc_fh_prs_10yr_quint_f ///
	crc_fh_10yr_quintile_f crc_fh_10yr_quint_f ///
	crc_prs_10yr_quintile_f crc_prs_10yr_quint_f ///
	crc_newfh_prs_10yr_quintile_f crc_newfh_prs_10yr_quint_f ///
	crc_newfh_10yr_quintile_f crc_newfh_10yr_quint_f ///
	crc_newprs_10yr_quintile_f crc_newprs_10yr_quint_f ///
	crc_new_10yr_quintile_f crc_new_10yr_quint_f ///
	crc_av_10yr_quintile_m crc_av_10yr_quint_m ///
	crc_fh_prs_10yr_quintile_m crc_fh_prs_10yr_quint_m ///
	crc_fh_10yr_quintile_m crc_fh_10yr_quint_m ///
	crc_prs_10yr_quintile_m crc_prs_10yr_quint_m ///
	crc_newfh_prs_10yr_quintile_m crc_newfh_prs_10yr_quint_m ///
	crc_newfh_10yr_quintile_m crc_newfh_10yr_quint_m ///
	crc_newprs_10yr_quintile_m crc_newprs_10yr_quint_m ///
	crc_new_10yr_quintile_m crc_new_10yr_quint_m 
	
*** deciles of risk ***

* women

pctile crc_av_10yr_decile_f = crc_av_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_av_10yr_dec_f = crc_av_10yr, ///
	cutpoints(crc_av_10yr_decile_f)
generate crc_av_10yr_dec=crc_av_10yr_dec_f if sex==0
label variable crc_av_10yr_dec "Average 10-year deciles"

pctile crc_fh_prs_10yr_decile_f = crc_fh_prs_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_fh_prs_10yr_dec_f = crc_fh_prs_10yr, ///
	cutpoints(crc_fh_prs_10yr_decile_f)
generate crc_fh_prs_10yr_dec=crc_fh_prs_10yr_dec_f if sex==0
label variable crc_fh_prs_10yr_dec "FH & PRS model 10-year deciles"

pctile crc_fh_10yr_decile_f = crc_fh_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_fh_10yr_dec_f = crc_fh_10yr, ///
	cutpoints(crc_fh_10yr_decile_f)
generate crc_fh_10yr_dec=crc_fh_10yr_dec_f if sex==0
label variable crc_fh_10yr_dec "FH model 10-year deciles"

pctile crc_prs_10yr_decile_f = crc_prs_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_prs_10yr_dec_f = crc_prs_10yr, ///
	cutpoints(crc_prs_10yr_decile_f)
generate crc_prs_10yr_dec=crc_prs_10yr_dec_f if sex==0
label variable crc_prs_10yr_dec "PRS model 10-year deciles"

pctile crc_newfh_prs_10yr_decile_f = crc_newfh_prs_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_newfh_prs_10yr_dec_f = crc_newfh_prs_10yr, ///
	cutpoints(crc_newfh_prs_10yr_decile_f)
generate crc_newfh_prs_10yr_dec=crc_newfh_prs_10yr_dec_f if sex==0
label variable crc_newfh_prs_10yr_dec "New FH & PRS model 10-year deciles"

pctile crc_newfh_10yr_decile_f = crc_newfh_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_newfh_10yr_dec_f = crc_newfh_10yr, cutpoints(crc_newfh_10yr_decile_f)
generate crc_newfh_10yr_dec=crc_newfh_10yr_dec_f if sex==0
label variable crc_newfh_10yr_dec "New FH model 10-year deciles"

pctile crc_newprs_10yr_decile_f = crc_newprs_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_newprs_10yr_dec_f = crc_newprs_10yr, ///
	cutpoints(crc_newprs_10yr_decile_f)
generate crc_newprs_10yr_dec=crc_newprs_10yr_dec_f if sex==0
label variable crc_newprs_10yr_dec "New PRS model 10-year deciles"

pctile crc_new_10yr_decile_f = crc_new_10yr if sex==0 & ///
	ancestry_code==9, nq(10)
xtile crc_new_10yr_dec_f = crc_new_10yr, cutpoints(crc_new_10yr_decile_f)
generate crc_new_10yr_dec=crc_new_10yr_dec_f if sex==0
label variable crc_new_10yr_dec "New model 10-year deciles"

* men
pctile crc_av_10yr_decile_m = crc_av_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_av_10yr_dec_m = crc_av_10yr, ///
	cutpoints(crc_av_10yr_decile_m)
replace crc_av_10yr_dec=crc_av_10yr_dec_m if sex==1

pctile crc_fh_prs_10yr_decile_m = crc_fh_prs_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_fh_prs_10yr_dec_m = crc_fh_prs_10yr, ///
	cutpoints(crc_fh_prs_10yr_decile_m)
replace crc_fh_prs_10yr_dec=crc_fh_prs_10yr_dec_m if sex==1

pctile crc_fh_10yr_decile_m = crc_fh_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_fh_10yr_dec_m = crc_fh_10yr, ///
	cutpoints(crc_fh_10yr_decile_m)
replace crc_fh_10yr_dec=crc_fh_10yr_dec_m if sex==1

pctile crc_prs_10yr_decile_m = crc_prs_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_prs_10yr_dec_m = crc_prs_10yr, ///
	cutpoints(crc_prs_10yr_decile_m)
replace crc_prs_10yr_dec=crc_prs_10yr_dec_m if sex==1

pctile crc_newfh_prs_10yr_decile_m = crc_newfh_prs_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_newfh_prs_10yr_dec_m = crc_newfh_prs_10yr, ///
	cutpoints(crc_newfh_prs_10yr_decile_m)
replace crc_newfh_prs_10yr_dec=crc_newfh_prs_10yr_dec_m if sex==1

pctile crc_newfh_10yr_decile_m = crc_newfh_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_newfh_10yr_dec_m = crc_newfh_10yr, cutpoints(crc_newfh_10yr_decile_m)
replace crc_newfh_10yr_dec=crc_newfh_10yr_dec_m if sex==1

pctile crc_newprs_10yr_decile_m = crc_newprs_10yr if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_newprs_10yr_dec_m = crc_newprs_10yr, ///
	cutpoints(crc_newprs_10yr_decile_m)
replace crc_newprs_10yr_dec=crc_newprs_10yr_dec_m if sex==1

pctile crc_new_10yr_decile_m = crc_new_10yr  if sex==1 & ///
	ancestry_code==9, nq(10)
xtile crc_new_10yr_dec_m = crc_new_10yr, cutpoints(crc_new_10yr_decile_m)
replace crc_new_10yr_dec=crc_new_10yr_dec_m if sex==1

drop crc_av_10yr_decile_f crc_av_10yr_dec_f ///
	crc_fh_prs_10yr_decile_f crc_fh_prs_10yr_dec_f ///
	crc_fh_10yr_decile_f crc_fh_10yr_dec_f ///
	crc_prs_10yr_decile_f crc_prs_10yr_dec_f ///
	crc_newfh_prs_10yr_decile_f crc_newfh_prs_10yr_dec_f ///
	crc_newfh_10yr_decile_f crc_newfh_10yr_dec_f ///
	crc_newprs_10yr_decile_f crc_newprs_10yr_dec_f ///
	crc_new_10yr_decile_f crc_new_10yr_dec_f ///
	crc_av_10yr_decile_m crc_av_10yr_dec_m ///
	crc_fh_prs_10yr_decile_m crc_fh_prs_10yr_dec_m ///
	crc_fh_10yr_decile_m crc_fh_10yr_dec_m ///
	crc_prs_10yr_decile_m crc_prs_10yr_dec_m ///
	crc_newfh_prs_10yr_decile_m crc_newfh_prs_10yr_dec_m ///
	crc_newfh_10yr_decile_m crc_newfh_10yr_dec_m ///
	crc_newprs_10yr_decile_m crc_newprs_10yr_dec_m ///
	crc_new_10yr_decile_m crc_new_10yr_dec_m 	
	

*** drop variables used to calculate risk

drop incid_a0-incid_s_life84 crc_av_10_a0-crc_av_10_a9 ///
	crc_fh_prs_10_a0-crc_fh_prs_10_a9 crc_fh_10_a0-crc_fh_10_a9 ///
	crc_prs_10_a0-crc_prs_10_a9 crc_newfh_prs_10_a0-crc_newfh_prs_10_a9 ///
	crc_newfh_10_a0-crc_newfh_10_a9 crc_newprs_10_a0-crc_newprs_10_a9 ///
	crc_new_10_a0-crc_new_10_a9 crc_av_life0-crc_av_life84 ///
	crc_fh_prs_life0-crc_fh_prs_life84 crc_fh_life0-crc_fh_life84 ///
	crc_prs_life0-crc_prs_life84 crc_newfh_prs_life0-crc_newfh_prs_life84 ///
	crc_newfh_life0-crc_newfh_life84 crc_newprs_life0-crc_newprs_life84 ///
	crc_new_life0-crc_new_life84

/*
*** 140-SNP PRS
colorpalette tab Hue circle	
twoway histogram prs_xb if ancestry_code==1, recast(line) lcolor("`r(p9)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==2, recast(line) lcolor("`r(p11)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==3, recast(line) lcolor("`r(p17)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==4, recast(line) lcolor("`r(p1)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==5, recast(line) lcolor("`r(p6)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==6, recast(line) lcolor("`r(p15)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==7, recast(line) lcolor("`r(p12)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==8, recast(line) lcolor("`r(p13)'") ///
		lwidth(medthick) || ///
	histogram prs_xb if ancestry_code==9, recast(line) lcolor(black) ///
		lwidth(medthick) ///
	xline(0, lcolor(gs8)) ///
	xscale(range(-2.25 2.25)) yscale(range(0 1.1)) ///
	xtitle("140-SNP PRS (lnRR)", height(5)) xlabel(-2(0.5)2, format(%3.1f)) ///
	ytitle("Density", height(5)) ylabel(0(0.2)1.0, angle(horizontal) ///
		format(%3.1f)) ///
	legend(region(lwidth(none)) size(*0.8) order(9 1 5 6 8 4 3 2 7) cols(3) ///
	label(1 "Ashkenazi") label(2 "Caribbean") label(3 "East Asia") ///
	label(4 "South Asia") label(5 "Middle East") label(6 "Southern Europe") ///
	label(7 "African") label(8 "Eastern Europe") label(9 "United Kingdom")) ///
	name(ancestry_140SNP, replace) 

graph save ancestry_140SNP "Graphs\PRS_distribution_xb.gph", replace
*graph export ethnicity_140SNP.tif, name(ancestry_140SNP) replace 

* summary statistics for the PRS distributions
tabstat	prs_rr, by(ancestry) statistics(N mean sd p50 iqr min max) ///
	columns(statistics)
*/
	
	
save "CR_testing_analysis.dta", replace
