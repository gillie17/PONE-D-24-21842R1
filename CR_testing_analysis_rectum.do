*** CRC model testing analysis in people of UK ancestry: rectum ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Testing_analysis_rectum_`filedate'.log", ///
	name(testing_rectum) replace

use "CR_testing_analysis.dta", clear

keep if ancestry_code==9 & cr_c10!=1 // exclude colon cases

* restricted to 10 years of follow-up
stset cr_r_cox_censor10_age, id(n_eid) enter(time age_calc) failure(cr_r10==1)


/* ANALYSIS OF UK ANCESTRY PEOPLE */

*** association per SD

* 10-year risk

stcox crc_av_10yr_logit_sd if sex==0
stcox crc_fh_10yr_logit_sd if sex==0
stcox crc_fh_prs_10yr_logit_sd if sex==0
stcox crc_newfh_prs_10yr_logit_sd if sex==0
stcox crc_new_10yr_logit_sd if sex==0

stcox crc_av_10yr_logit_sd if sex==1
stcox crc_fh_10yr_logit_sd if sex==1
stcox crc_fh_prs_10yr_logit_sd if sex==1
stcox crc_newfh_prs_10yr_logit_sd if sex==1
stcox crc_new_10yr_logit_sd if sex==1


*** discrimination

* 10-year risk

* women
somersd cr_r10 crc_av_10yr crc_fh_10yr crc_fh_prs_10yr crc_newfh_prs_10yr ///
	crc_new_10yr if sex==0, transf(c)

test crc_av_10yr=0.5 
test crc_fh_10yr=0.5 
test crc_fh_prs_10yr=0.5 
test crc_newfh_prs_10yr=0.5 
test crc_new_10yr=0.5 

lincom crc_new_10yr-crc_fh_prs_10yr 
lincom crc_new_10yr-crc_newfh_prs_10yr
lincom crc_newfh_prs_10yr-crc_fh_prs_10yr

lincom crc_fh_prs_10yr-crc_av_10yr 
lincom crc_newfh_prs_10yr-crc_av_10yr 
lincom crc_new_10yr-crc_av_10yr 

lincom crc_fh_prs_10yr-crc_fh_10yr 
lincom crc_newfh_prs_10yr-crc_fh_10yr 
lincom crc_new_10yr-crc_fh_10yr 


* men
somersd cr_r10 crc_av_10yr crc_fh_10yr crc_fh_prs_10yr crc_newfh_prs_10yr ///
	crc_new_10yr if sex==1, transf(c)

test crc_av_10yr=0.5 
test crc_fh_10yr=0.5 
test crc_fh_prs_10yr=0.5 
test crc_newfh_prs_10yr=0.5 
test crc_new_10yr=0.5 


lincom crc_new_10yr-crc_fh_prs_10yr 
lincom crc_new_10yr-crc_newfh_prs_10yr
lincom crc_newfh_prs_10yr-crc_fh_prs_10yr

lincom crc_fh_prs_10yr-crc_av_10yr 
lincom crc_newfh_prs_10yr-crc_av_10yr 
lincom crc_new_10yr-crc_av_10yr 

lincom crc_fh_prs_10yr-crc_fh_10yr 
lincom crc_newfh_prs_10yr-crc_fh_10yr 
lincom crc_new_10yr-crc_fh_10yr 


*** calibration 

* slope

* women
logit cr_r10 crc_av_10yr_logit if sex==0
test crc_av_10yr_logit=1

logit cr_r10 crc_fh_10yr_logit if sex==0
test crc_fh_10yr_logit=1

logit cr_r10 crc_prs_10yr_logit if sex==0
test crc_prs_10yr_logit=1

logit cr_r10 crc_fh_prs_10yr_logit if sex==0
test crc_fh_prs_10yr_logit=1
estimates store fhprsf_slope

logit cr_r10 crc_newfh_10yr_logit if sex==0
test crc_newfh_10yr_logit=1

logit cr_r10 crc_newprs_10yr_logit if sex==0
test crc_newprs_10yr_logit=1

logit cr_r10 crc_newfh_prs_10yr_logit if sex==0
test crc_newfh_prs_10yr_logit=1
estimates store newfhprsf_slope

logit cr_r10 crc_new_10yr_logit if sex==0
test crc_new_10yr_logit=1
estimates store newf_slope

suest fhprsf_slope newf_slope
test [fhprsf_slope_cr_r10]crc_fh_prs_10yr_logit = ///
	[newf_slope_cr_r10]crc_new_10yr_logit 

suest newfhprsf_slope newf_slope
test [newfhprsf_slope_cr_r10]crc_newfh_prs_10yr_logit = ///
	[newf_slope_cr_r10]crc_new_10yr_logit 

* men
logit cr_r10 crc_av_10yr_logit if sex==1
test crc_av_10yr_logit=1

logit cr_r10 crc_fh_10yr_logit if sex==1
test crc_fh_10yr_logit=1

logit cr_r10 crc_prs_10yr_logit if sex==1
test crc_prs_10yr_logit=1

logit cr_r10 crc_fh_prs_10yr_logit if sex==1
test crc_fh_prs_10yr_logit=1
estimates store fhprsm_slope

logit cr_r10 crc_newfh_10yr_logit if sex==1
test crc_newfh_10yr_logit=1

logit cr_r10 crc_newprs_10yr_logit if sex==1
test crc_newprs_10yr_logit=1

logit cr_r10 crc_newfh_prs_10yr_logit if sex==1
test crc_newfh_prs_10yr_logit=1
estimates store newfhprsm_slope

logit cr_r10 crc_new_10yr_logit if sex==1
test crc_new_10yr_logit=1
estimates store newm_slope

suest fhprsm_slope newm_slope
test [fhprsm_slope_cr_r10]crc_fh_prs_10yr_logit = ///
	[newm_slope_cr_r10]crc_new_10yr_logit 

suest newfhprsm_slope newm_slope
test [newfhprsm_slope_cr_r10]crc_newfh_prs_10yr_logit = ///
	[newm_slope_cr_r10]crc_new_10yr_logit 

* intercept

* women
constraint define 1 crc_av_10yr_logit=1
logit cr_r10 crc_av_10yr_logit if sex==0, constraints(1) 
estimates store avf_intercept

constraint define 2 crc_fh_10yr_logit=1
logit cr_r10 crc_fh_10yr_logit if sex==0, constraints(2) 
estimates store fhf_intercept

constraint define 3 crc_prs_10yr_logit=1
logit cr_r10 crc_prs_10yr_logit if sex==0, constraints(3) 

constraint define 4 crc_fh_prs_10yr_logit=1
logit cr_r10 crc_fh_prs_10yr_logit if sex==0, constraints(4) 
estimates store fhprsf_intercept

constraint define 5 crc_newfh_10yr_logit=1
logit cr_r10 crc_newfh_10yr_logit if sex==0, constraints(5) 

constraint define 6 crc_newprs_10yr_logit=1
logit cr_r10 crc_newprs_10yr_logit if sex==0, constraints(6) 

constraint define 7 crc_newfh_prs_10yr_logit=1
logit cr_r10 crc_newfh_prs_10yr_logit if sex==0, constraints(7) 
estimates store newfhprsf_intercept

constraint define 8 crc_new_10yr_logit=1
logit cr_r10 crc_new_10yr_logit if sex==0, constraints(8) 
estimates store newf_intercept

suest newfhprsf_intercept newf_intercept
test [newfhprsf_intercept_cr_r10]_cons = [newf_intercept_cr_r10]_cons 

suest avf_intercept fhf_intercept
test [avf_intercept_cr_r10]_cons = [fhf_intercept_cr_r10]_cons 

suest avf_intercept fhprsf_intercept
test [avf_intercept_cr_r10]_cons = [fhprsf_intercept_cr_r10]_cons 

suest avf_intercept newfhprsf_intercept
test [avf_intercept_cr_r10]_cons = [newfhprsf_intercept_cr_r10]_cons 

suest avf_intercept newf_intercept
test [avf_intercept_cr_r10]_cons = [newf_intercept_cr_r10]_cons 


* men

logit cr_r10 crc_av_10yr_logit if sex==1, constraints(1) 
estimates store avm_intercept

logit cr_r10 crc_fh_10yr_logit if sex==1, constraints(2) 
estimates store fhm_intercept

logit cr_r10 crc_prs_10yr_logit if sex==1, constraints(3) 

logit cr_r10 crc_fh_prs_10yr_logit if sex==1, constraints(4) 
estimates store fhprsm_intercept

logit cr_r10 crc_newfh_10yr_logit if sex==1, constraints(5) 

logit cr_r10 crc_newprs_10yr_logit if sex==1, constraints(6) 

logit cr_r10 crc_newfh_prs_10yr_logit if sex==1, constraints(7) 
estimates store newfhprsm_intercept

logit cr_r10 crc_new_10yr_logit if sex==1, constraints(8) 
estimates store newm_intercept


suest newfhprsm_intercept newm_intercept
test [newfhprsm_intercept_cr_r10]_cons = [newm_intercept_cr_r10]_cons 

suest avm_intercept fhm_intercept
test [avm_intercept_cr_r10]_cons = [fhm_intercept_cr_r10]_cons 

suest avm_intercept fhprsm_intercept
test [avm_intercept_cr_r10]_cons = [fhprsm_intercept_cr_r10]_cons 

suest avm_intercept newfhprsm_intercept
test [avm_intercept_cr_r10]_cons = [newfhprsm_intercept_cr_r10]_cons 

suest avm_intercept newm_intercept
test [avm_intercept_cr_r10]_cons = [newm_intercept_cr_r10]_cons 


log close testing_rectum