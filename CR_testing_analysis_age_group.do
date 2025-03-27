*** CRC model testing analysis in people of UK ancestry ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Testing_analysis_agegp_`filedate'.log", ///
	name(testing_agegp) replace

use "CR_testing_analysis.dta", clear

keep if ancestry_code==9

* restricted to 10 years of follow-up
stset cr_cox_censor10_age, id(n_eid) enter(time age_calc) failure(cr10==1)


/* ANALYSIS OF UK ANCESTRY PEOPLE */

	
**# association per SD

foreach n in 0 1 {

di ""
di "Sex = " `n'
di ""

foreach m in 1 2 3 {

di ""
di " age group = `m'"

stcox crc_av_10yr_logit_sd if sex==`n' & age_gp10==`m'
stcox crc_fh_10yr_logit_sd if sex==`n' & age_gp10==`m'
stcox crc_fh_prs_10yr_logit_sd if sex==`n' & age_gp10==`m'
stcox crc_newfh_prs_10yr_logit_sd if sex==`n' & age_gp10==`m'
stcox crc_new_10yr_logit_sd if sex==`n' & age_gp10==`m'
}
}

* women vs men

foreach m in 1 2 3 {

di ""
di " age group = `m'"

stcox c.crc_av_10yr_logit_sd##i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_av_10yr_logit_sd] = _b[1.sex#c.crc_av_10yr_logit_sd]

stcox c.crc_fh_10yr_logit_sd##i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_fh_10yr_logit_sd] = _b[1.sex#c.crc_fh_10yr_logit_sd]

stcox c.crc_fh_prs_10yr_logit_sd##i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_fh_prs_10yr_logit_sd] = _b[1.sex#c.crc_fh_prs_10yr_logit_sd]

stcox c.crc_newfh_prs_10yr_logit_sd##i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_newfh_prs_10yr_logit_sd] = _b[1.sex#c.crc_newfh_prs_10yr_logit_sd]

stcox c.crc_new_10yr_logit_sd##i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_new_10yr_logit_sd] = _b[1.sex#c.crc_new_10yr_logit_sd]

}
	
	
*** test difference between HRs for different models

* do poisson version of cox then use suest

foreach n in 0 1 {
foreach m in 1 2 3 {

preserve

keep if sex==`n'
keep if age_gp10==`m'

* split by age (faster than splitting at failures)

stsplit agegp_1yr, at(40(1)80) 
generate time_exp = _t - _t0

foreach s in av fh fh_prs newfh_prs new {

poisson _d ibn.agegp_1yr crc_`s'_10yr_logit_sd, exposure(time_exp) noconstant 
estimates store `s'
}

di ""
di "Sex = " `n'
di ""
di " age group = `m'"

quietly suest av newfh_prs
test [av__d]crc_av_10yr_logit_sd = [newfh_prs__d]crc_newfh_prs_10yr_logit_sd 

quietly suest fh newfh_prs
test [fh__d]crc_fh_10yr_logit_sd = [newfh_prs__d]crc_newfh_prs_10yr_logit_sd 

quietly suest fh_prs newfh_prs
test [fh_prs__d]crc_fh_prs_10yr_logit_sd = [newfh_prs__d]crc_newfh_prs_10yr_logit_sd 

quietly suest av new
test [av__d]crc_av_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

quietly suest fh new
test [fh__d]crc_fh_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

quietly suest fh_prs new
test [fh_prs__d]crc_fh_prs_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

quietly suest newfh_prs new
test [newfh_prs__d]crc_newfh_prs_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

restore

}
}


*** test difference between HRs for different age groups

foreach n in 0 1 {
foreach s in av fh fh_prs newfh_prs new  {

preserve

keep if sex==`n'

* split by age (faster than splitting at failures)

stsplit agegp_1yr, at(40(1)80) 
generate time_exp = _t - _t0

foreach m in 1 2 3 {

poisson _d ibn.agegp_1yr crc_`s'_10yr_logit_sd if age_gp10==`m', exposure(time_exp) noconstant 
estimates store `s'_`m'
}

quietly suest `s'_1 `s'_2
test [`s'_1__d]crc_`s'_10yr_logit_sd = [`s'_2__d]crc_`s'_10yr_logit_sd 

quietly suest `s'_2 `s'_3
test [`s'_2__d]crc_`s'_10yr_logit_sd = [`s'_3__d]crc_`s'_10yr_logit_sd 

quietly suest `s'_1 `s'_3
test [`s'_1__d]crc_`s'_10yr_logit_sd = [`s'_3__d]crc_`s'_10yr_logit_sd 
di " "

restore

}
}


**# calibration 

*** slope

foreach n in 0 1 {
foreach m in 1 2 3 {

di ""
di "Sex = " `n'
di ""
di "age group = `m'"


logit cr10 crc_av_10yr_logit if sex==`n' & age_gp10==`m'
test crc_av_10yr_logit=1
estimates store av_slope

logit cr10 crc_fh_10yr_logit if sex==`n' & age_gp10==`m'
test crc_fh_10yr_logit=1
estimates store fh_slope

logit cr10 crc_fh_prs_10yr_logit if sex==`n' & age_gp10==`m'
test crc_fh_prs_10yr_logit=1
estimates store fh_prs_slope

logit cr10 crc_newfh_prs_10yr_logit if sex==`n' & age_gp10==`m'
test crc_newfh_prs_10yr_logit=1
estimates store newfh_prs_slope

logit cr10 crc_new_10yr_logit if sex==`n' & age_gp10==`m'
test crc_new_10yr_logit=1
estimates store new_slope


quietly suest av_slope newfh_prs_slope
test [av_slope_cr10]crc_av_10yr_logit = [newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit 

quietly suest fh_slope newfh_prs_slope
test [fh_slope_cr10]crc_fh_10yr_logit = [newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit 

quietly suest fh_prs_slope newfh_prs_slope
test [fh_prs_slope_cr10]crc_fh_prs_10yr_logit = [newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit 

quietly suest av_slope new_slope
test [av_slope_cr10]crc_av_10yr_logit = [new_slope_cr10]crc_new_10yr_logit 

quietly suest fh_slope new_slope
test [fh_slope_cr10]crc_fh_10yr_logit = [new_slope_cr10]crc_new_10yr_logit 

quietly suest fh_prs_slope new_slope
test [fh_prs_slope_cr10]crc_fh_prs_10yr_logit = [new_slope_cr10]crc_new_10yr_logit 

quietly suest newfh_prs_slope new_slope
test [newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit = [new_slope_cr10]crc_new_10yr_logit 
}
}

* test sex interaction

foreach m in 1 2 3 {

logit cr10 c.crc_av_10yr_logit#i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_av_10yr_logit] = _b[1.sex#c.crc_av_10yr_logit]

logit cr10 c.crc_fh_10yr_logit#i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_fh_10yr_logit] = _b[1.sex#c.crc_fh_10yr_logit]

logit cr10 c.crc_fh_prs_10yr_logit#i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_fh_prs_10yr_logit] = _b[1.sex#c.crc_fh_prs_10yr_logit]
	
logit cr10 c.crc_newfh_prs_10yr_logit#i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_newfh_prs_10yr_logit] = _b[1.sex#c.crc_newfh_prs_10yr_logit]

logit cr10 c.crc_new_10yr_logit#i.sex if age_gp10==`m'
test _b[0b.sex#c.crc_new_10yr_logit] = _b[1.sex#c.crc_new_10yr_logit]

}


*** test difference between HRs for different age groups

foreach n in 0 1 {
foreach s in av fh fh_prs newfh_prs new  {

preserve

keep if sex==`n'

logit cr10 crc_`s'_10yr_logit if age_gp10==1
estimates store `s'_slope_1

logit cr10 crc_`s'_10yr_logit if age_gp10==2
estimates store `s'_slope_2

logit cr10 crc_`s'_10yr_logit if age_gp10==3
estimates store `s'_slope_3

suest `s'_slope_1 `s'_slope_2
test [`s'_slope_1_cr10]crc_`s'_10yr_logit = [`s'_slope_2_cr10]crc_`s'_10yr_logit

suest `s'_slope_2 `s'_slope_3
test [`s'_slope_2_cr10]crc_`s'_10yr_logit = [`s'_slope_3_cr10]crc_`s'_10yr_logit

suest `s'_slope_1 `s'_slope_3
test [`s'_slope_1_cr10]crc_`s'_10yr_logit = [`s'_slope_3_cr10]crc_`s'_10yr_logit

restore

}
}



*** intercept

foreach m in 1 2 3 {
foreach n in 0 1 {

di ""
di "Sex =  `n'"
di ""
di "age group = `m'"

constraint define 1 crc_av_10yr_logit=1
logit cr10 crc_av_10yr_logit if sex==`n' & age_gp10==`m', constraints(1) 
estimates store av_intercept_`n'

constraint define 2 crc_fh_10yr_logit=1
logit cr10 crc_fh_10yr_logit if sex==`n' & age_gp10==`m', constraints(2) 
estimates store fh_intercept_`n'

constraint define 3 crc_fh_prs_10yr_logit=1
logit cr10 crc_fh_prs_10yr_logit if sex==`n' & age_gp10==`m', constraints(3) 
estimates store fh_prs_intercept_`n'

constraint define 4 crc_newfh_prs_10yr_logit=1
logit cr10 crc_newfh_prs_10yr_logit if sex==`n' & age_gp10==`m', constraints(4) 
estimates store newfh_prs_intercept_`n'

constraint define 5 crc_new_10yr_logit=1
logit cr10 crc_new_10yr_logit if sex==`n' & age_gp10==`m', constraints(5) 
estimates store new_intercept_`n'

suest av_intercept_`n' newfh_prs_intercept_`n'
test [av_intercept_`n'_cr10]_cons = [newfh_prs_intercept_`n'_cr10]_cons 

suest fh_intercept_`n' newfh_prs_intercept_`n'
test [fh_intercept_`n'_cr10]_cons = [newfh_prs_intercept_`n'_cr10]_cons 

suest fh_prs_intercept_`n' newfh_prs_intercept_`n'
test [fh_prs_intercept_`n'_cr10]_cons = [newfh_prs_intercept_`n'_cr10]_cons 

suest av_intercept_`n' new_intercept_`n'
test [av_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest fh_intercept_`n' new_intercept_`n'
test [fh_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest fh_prs_intercept_`n' new_intercept_`n'
test [fh_prs_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest newfh_prs_intercept_`n' new_intercept_`n'
test [newfh_prs_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

}

* test sex interaction

suest av_intercept_0 av_intercept_1
test [av_intercept_0_cr10]_cons = [av_intercept_1_cr10]_cons 

suest fh_intercept_0 fh_intercept_1
test [fh_intercept_0_cr10]_cons = [fh_intercept_1_cr10]_cons 

suest fh_prs_intercept_0 fh_prs_intercept_1
test [fh_prs_intercept_0_cr10]_cons = [fh_prs_intercept_1_cr10]_cons 

suest newfh_prs_intercept_0 newfh_prs_intercept_1
test [newfh_prs_intercept_0_cr10]_cons = [newfh_prs_intercept_1_cr10]_cons 

suest new_intercept_0 new_intercept_1
test [new_intercept_0_cr10]_cons = [new_intercept_1_cr10]_cons 

}

	
*** test difference between intercepts for different age groups

foreach n in 0 1 {
foreach s in av fh fh_prs newfh_prs new  {

preserve

keep if sex==`n'

constraint define 1 crc_`s'_10yr_logit=1

logit cr10 crc_`s'_10yr_logit if age_gp10==1, constraints(1) 
estimates store `s'_intercept_1

logit cr10 crc_`s'_10yr_logit if age_gp10==2, constraints(1) 
estimates store `s'_intercept_2

logit cr10 crc_`s'_10yr_logit if age_gp10==3, constraints(1) 
estimates store `s'_intercept_3


suest `s'_intercept_1 `s'_intercept_2
test [`s'_intercept_1_cr10]_cons = [`s'_intercept_2_cr10]_cons

suest `s'_intercept_2 `s'_intercept_3
test [`s'_intercept_2_cr10]_cons = [`s'_intercept_3_cr10]_cons

suest `s'_intercept_1 `s'_intercept_3
test [`s'_intercept_1_cr10]_cons = [`s'_intercept_3_cr10]_cons



restore

}
}
	
	
	
log close testing_agegp
exit