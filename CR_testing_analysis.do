*** CRC model testing analysis in people of UK ancestry ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Testing_analysis_`filedate'.log", ///
	name(testing) replace

use "CR_testing_analysis.dta", clear

keep if ancestry_code==9

* restricted to 10 years of follow-up
stset cr_cox_censor10_age, id(n_eid) enter(time age_calc) failure(cr10==1)


/* ANALYSIS OF UK ANCESTRY PEOPLE */


**# summary stats ***

tabstat crc_fh_rr crc_fh_prs_rr crc_newfh_prs_rr crc_new_rr, by(sex) ///
	statistics(mean sd p50 iqr min max) ///
	columns(statistics) nototal

generate crc_fh_xb=ln(crc_fh_rr)
generate crc_fh_prs_xb=ln(crc_fh_prs_rr)
generate crc_newfh_prs_xb=ln(crc_newfh_prs_rr)
generate crc_new_xb=ln(crc_new_rr)
	
tabstat crc_fh_xb crc_fh_prs_xb crc_newfh_prs_xb crc_new_xb, by(sex) ///
	statistics(mean sd p50 iqr min max) ///
	columns(statistics) nototal	
	
tabstat crc_av_10yr crc_fh_10yr crc_fh_prs_10yr crc_newfh_prs_10yr ///
	crc_new_10yr, by(sex) statistics(mean sd p50 iqr min max) ///
	columns(statistics) nototal

	
**# association per SD

foreach n in 0 1 {

di ""
di "Sex = " `n'
di ""

foreach v of varlist crc_av_10yr_logit_sd crc_fh_10yr_logit_sd ///
	crc_fh_prs_10yr_logit_sd crc_newfh_prs_10yr_logit_sd crc_new_10yr_logit_sd {

stcox `v' if sex==`n'

}
}

* women vs men

foreach v of varlist crc_av_10yr_logit_sd crc_fh_10yr_logit_sd ///
	crc_fh_prs_10yr_logit_sd crc_newfh_prs_10yr_logit_sd crc_new_10yr_logit_sd {

stcox c.`v'##i.sex
test _b[0b.sex#c.`v'] = _b[1.sex#c.`v']

}
	
	
*** test difference between HRs for different models

* do poisson version of cox then use suest

foreach n in 0 1 {

preserve

keep if sex==`n'

di ""
di "Sex = " `n'
di ""

* split by age (faster than splitting at failures)

stsplit agegp_1yr, at(40(1)80) 
generate time_exp = _t - _t0

foreach s in av fh fh_prs newfh_prs new {

poisson _d ibn.agegp_1yr crc_`s'_10yr_logit_sd, exposure(time_exp) noconstant 
estimates store `s'
}

quietly suest av fh_prs
test [av__d]crc_av_10yr_logit_sd = [fh_prs__d]crc_fh_prs_10yr_logit_sd 

quietly suest av newfh_prs
test [av__d]crc_av_10yr_logit_sd = [newfh_prs__d]crc_newfh_prs_10yr_logit_sd 

quietly suest av new
test [av__d]crc_av_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

quietly suest fh fh_prs
test [fh__d]crc_fh_10yr_logit_sd = [fh_prs__d]crc_fh_prs_10yr_logit_sd 

quietly suest fh newfh_prs
test [fh__d]crc_fh_10yr_logit_sd = [newfh_prs__d]crc_newfh_prs_10yr_logit_sd 

quietly suest fh new
test [fh__d]crc_fh_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 


quietly suest newfh_prs new
test [newfh_prs__d]crc_newfh_prs_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

quietly suest fh_prs new
test [fh_prs__d]crc_fh_prs_10yr_logit_sd = [new__d]crc_new_10yr_logit_sd 

quietly suest fh_prs newfh_prs
test [fh_prs__d]crc_fh_prs_10yr_logit_sd = [newfh_prs__d]crc_newfh_prs_10yr_logit_sd 


restore

}


**# discrimination

foreach n in 0 1 {

di ""
di "Sex = " `n'
di "

somersd cr10 crc_av_10yr crc_fh_10yr crc_fh_prs_10yr crc_newfh_prs_10yr ///
	crc_new_10yr if sex==`n', transf(c)
	
estimates store harrell_`n'	

test crc_av_10yr=0.5 
test crc_fh_10yr=0.5 
test crc_fh_prs_10yr=0.5 
test crc_newfh_prs_10yr=0.5 
test crc_new_10yr=0.5 

lincom crc_fh_prs_10yr-crc_av_10yr 
lincom crc_fh_prs_10yr-crc_fh_10yr 

lincom crc_newfh_prs_10yr-crc_av_10yr 
lincom crc_newfh_prs_10yr-crc_fh_10yr

lincom crc_new_10yr-crc_av_10yr 
lincom crc_new_10yr-crc_fh_10yr

lincom crc_new_10yr-crc_fh_prs_10yr 
lincom crc_new_10yr-crc_newfh_prs_10yr
lincom crc_newfh_prs_10yr-crc_fh_prs_10yr

}

* z-test for comparisons between sex

foreach v of varlist crc_av_10yr crc_fh_10yr crc_fh_prs_10yr ///
	crc_newfh_prs_10yr crc_new_10yr {
foreach n in 0 1 {

di ""
di "Risk score = `v'"
di "Sex = " `n'
di ""

somersd cr10 `v' if sex==`n', transf(c)
scalar b_`v'_`n'=e(b)[1,1]
scalar v_`v'_`n'=e(V)[1,1]	
	
}
}

foreach v of varlist crc_av_10yr crc_fh_10yr crc_fh_prs_10yr ///
	crc_newfh_prs_10yr crc_new_10yr {

scalar z_`v'=abs(b_`v'_0 - b_`v'_1)/sqrt(v_`v'_0 + v_`v'_1)

di ""
di "`v' women vs men: " _column(35) "z = " z_`v' ///
	_column(55) "P = " (1 - normal(z_`v'))*2

di ""
}



**# calibration 

*** slope

foreach n in 0 1 {

di ""
di "Sex = " `n'
di "

logit cr10 crc_av_10yr_logit if sex==`n'
test crc_av_10yr_logit=1
estimates store av_slope

logit cr10 crc_fh_10yr_logit if sex==`n'
test crc_fh_10yr_logit=1
estimates store fh_slope

logit cr10 crc_fh_prs_10yr_logit if sex==`n'
test crc_fh_prs_10yr_logit=1
estimates store fh_prs_slope

logit cr10 crc_newfh_prs_10yr_logit if sex==`n'
test crc_newfh_prs_10yr_logit=1
estimates store newfh_prs_slope

logit cr10 crc_new_10yr_logit if sex==`n'
test crc_new_10yr_logit=1
estimates store new_slope

quietly suest av_slope new_slope
test [av_slope_cr10]crc_av_10yr_logit = ///
	[new_slope_cr10]crc_new_10yr_logit 

quietly suest fh_slope new_slope
test [fh_slope_cr10]crc_fh_10yr_logit = ///
	[new_slope_cr10]crc_new_10yr_logit 

quietly suest fh_prs_slope new_slope
test [fh_prs_slope_cr10]crc_fh_prs_10yr_logit = ///
	[new_slope_cr10]crc_new_10yr_logit 

quietly suest av_slope newfh_prs_slope
test [av_slope_cr10]crc_av_10yr_logit = ///
	[newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit 
	
quietly suest fh_slope newfh_prs_slope
test [fh_slope_cr10]crc_fh_10yr_logit = ///
	[newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit 

quietly suest fh_prs_slope newfh_prs_slope
test [fh_prs_slope_cr10]crc_fh_prs_10yr_logit = ///
	[newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit 
	
quietly suest newfh_prs_slope new_slope
test [newfh_prs_slope_cr10]crc_newfh_prs_10yr_logit = ///
	[new_slope_cr10]crc_new_10yr_logit 
}

* test sex interaction

foreach v of varlist crc_av_10yr crc_fh_10yr crc_fh_prs_10yr ///
	crc_newfh_prs_10yr crc_new_10yr {

logit cr10 c.`v'_logit#i.sex 
test  _b[0b.sex#c.`v'_logit] = _b[1.sex#c.`v'_logit]

}

*** intercept

foreach n in 0 1 {

di ""
di "Sex = " `n'
di "

constraint define 1 crc_av_10yr_logit=1
logit cr10 crc_av_10yr_logit if sex==`n', constraints(1) 
estimates store av_intercept_`n'

constraint define 2 crc_fh_10yr_logit=1
logit cr10 crc_fh_10yr_logit if sex==`n', constraints(2) 
estimates store fh_intercept_`n'

constraint define 3 crc_fh_prs_10yr_logit=1
logit cr10 crc_fh_prs_10yr_logit if sex==`n', constraints(3) 
estimates store fh_prs_intercept_`n'

constraint define 4 crc_newfh_prs_10yr_logit=1
logit cr10 crc_newfh_prs_10yr_logit if sex==`n', constraints(4) 
estimates store newfh_prs_intercept_`n'

constraint define 5 crc_new_10yr_logit=1
logit cr10 crc_new_10yr_logit if sex==`n', constraints(5) 
estimates store new_intercept_`n'

suest av_intercept_`n' new_intercept_`n'
test [av_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest fh_intercept_`n' new_intercept_`n'
test [fh_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest fh_prs_intercept_`n' new_intercept_`n'
test [fh_prs_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest newfh_prs_intercept_`n' new_intercept_`n'
test [newfh_prs_intercept_`n'_cr10]_cons = [new_intercept_`n'_cr10]_cons 

suest av_intercept_`n' newfh_prs_intercept_`n'
test [av_intercept_`n'_cr10]_cons = [newfh_prs_intercept_`n'_cr10]_cons 

suest fh_intercept_`n' newfh_prs_intercept_`n'
test [fh_intercept_`n'_cr10]_cons = [newfh_prs_intercept_`n'_cr10]_cons 

suest fh_prs_intercept_`n' newfh_prs_intercept_`n'
test [fh_prs_intercept_`n'_cr10]_cons = [newfh_prs_intercept_`n'_cr10]_cons 
 
}

* test sex interaction

foreach v in av_intercept fh_intercept fh_prs_intercept ///
	newfh_prs_intercept new_intercept {

suest `v'_0 `v'_1
test [`v'_0_cr10]_cons = [`v'_1_cr10]_cons 

}


**# calibration plots 
* (clplot_edit.grec edits the reference line and graph colours)

*** women

* Average
pmcalplot crc_av_10yr cr10 if sex==0, ci bin(10) range(0 0.05) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title(" " "Average risks: women", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_av_f, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:A}", size(*1.5) just(left))

graph save "cplot_av_f" "Graphs\calibplot_av_women.gph", replace
graph export "Graphs\calibplot_av_women.tif", as(tif) ///
	name("cplot_av_f") replace


* Aviv's FH 
pmcalplot crc_fh_10yr cr10 if sex==0, ci bin(10) range(0 0.05) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("Family history alone model:" "women", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_fh_f, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:B}", size(*1.5) just(left))

graph save "cplot_fh_f" "Graphs\calibplot_fh_women.gph", replace
graph export "Graphs\calibplot_fh_women.tif", as(tif) ///
	name("cplot_fh_f") replace

* Aviv's FH and PRS
pmcalplot crc_fh_prs_10yr cr10 if sex==0, ci bin(10) range(0 0.05) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("Current family history and" "PRS model: women", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_fh_prs_f, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:C}", size(*1.5) just(left))

graph save "cplot_fh_prs_f" "Graphs\calibplot_fh_prs_women.gph", replace
graph export "Graphs\calibplot_fh_prs_women.tif", as(tif) ///
	name("cplot_fh_prs_f") replace

* new FH and PRS
pmcalplot crc_newfh_prs_10yr cr10 if sex==0, ci bin(10) range(0 0.05) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("New family history and " "PRS model: women", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_newfh_prs_f, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:D}", size(*1.5) just(left))

graph save "cplot_newfh_prs_f" "Graphs\calibplot_newfh_prs_women.gph", replace
graph export "Graphs\calibplot_newfh_prs_women.tif", as(tif) ///
	name("cplot_newfh_prs_f") replace	
	
* new
pmcalplot crc_new_10yr cr10 if sex==0, ci bin(10) range(0 0.05) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("New multivariable model:" "women", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_new_f, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:E}", size(*1.5) just(left))

graph save "cplot_new_f" "Graphs\calibplot_new_women.gph", replace
graph export "Graphs\calibplot_new_women.tif", as(tif) ///
	name("cplot_new_f") replace

graph combine "Graphs\calibplot_av_women.gph" "Graphs\calibplot_fh_women.gph" ///
	"Graphs\calibplot_fh_prs_women.gph" "Graphs\calibplot_newfh_prs_women.gph" ///
	"Graphs\calibplot_new_women.gph" , cols(2) ysize(9) xsize(5.8) iscale(0.65) ///
	graphregion(margin(small)) name(calib_women, replace)

graph save calib_women "Graphs\Calibration_women_combined.gph", replace
graph export "Graphs\Calibration_women_combined.tif", name(calib_women) replace	
	
	
*** men

* Average 
pmcalplot crc_av_10yr cr10 if sex==1, ci bin(10) range(0 0.0525) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title(" " "Average risks: men", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_av_m, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:A}", size(*1.5) just(left))

graph save "cplot_av_m" "Graphs\calibplot_av_men.gph", replace
graph export "Graphs\calibplot_av_men.tif", as(tif) ///
	name("cplot_av_m") replace


* Aviv's FH 
pmcalplot crc_fh_10yr cr10 if sex==1, ci bin(10) range(0 0.0525) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("Family history alone model:" "men", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_fh_m, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:B}", size(*1.5) just(left))

graph save "cplot_fh_m" "Graphs\calibplot_fh_men.gph", replace
graph export "Graphs\calibplot_fh_men.tif", as(tif) ///
	name("cplot_fh_m") replace

* Aviv's FH and PRS
pmcalplot crc_fh_prs_10yr cr10 if sex==1, ci bin(10) range(0 0.0525) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("Current family history and" "PRS model: men", size(*0.8)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_fh_prs_m, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:C}", size(*1.5) just(left))

graph save "cplot_fh_prs_m" "Graphs\calibplot_fh_prs_men.gph", replace
graph export "Graphs\calibplot_fh_prs_men.tif", as(tif) ///
	name("cplot_fh_prs_m") replace


* new FH and PRS
pmcalplot crc_newfh_prs_10yr cr10 if sex==1, ci bin(10) range(0 0.0525) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title("New family history and" "PRS model: men", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_newfh_prs_m, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:D}", size(*1.5) just(left))

graph save "cplot_newfh_prs_m" "Graphs\calibplot_newfh_prs_men.gph", replace
graph export "Graphs\calibplot_newfh_prs_men.tif", as(tif) ///
	name("cplot_newfh_prs_m") replace
	
* new
pmcalplot crc_new_10yr cr10 if sex==1, ci bin(10) range(0 0.0525) ///
	xsize(4) nospike nolowess nostatistics ///
	ciopts(lcolor(purple)) scatteropts(msymbol(O) mcolor(purple)) ///
	title(" " "New multivariable model:" "men", size(*0.7)) ///
	xtitle("Expected 10-year risk (%)", height(5)) ///
	ytitle("Observed 10-year risk (%)", height(4)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5") ///
	ylabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5", angle(h)) ///
	legend(order(2 3 1) label(1 "y = x") label(2 "Decile of risk") ///
		label(3 "95% CI") region(lwidth(none)) cols(1) size(*0.8) ///
		symxsize(5) keygap(*1) rowgap(*0.2) ring(0) position(5)) ///
		name(cplot_new_m, replace) play(clplot_edit) ///
	text(0.0575 -0.01 "{bf:E}", size(*1.5) just(left))

graph save "cplot_new_m" "Graphs\calibplot_new_men.gph", replace
graph export "Graphs\calibplot_new_men.tif", as(tif) ///
	name("cplot_new_m") replace

	
graph combine "Graphs\calibplot_av_men.gph" "Graphs\calibplot_fh_men.gph" ///
	"Graphs\calibplot_fh_prs_men.gph" "Graphs\calibplot_newfh_prs_men.gph" ///
	"Graphs\calibplot_new_men.gph" , cols(2) ysize(9) xsize(5.8) iscale(0.6) ///
	graphregion(margin(small)) name(calib_men, replace)

graph save calib_men "Graphs\Calibration_men_combined.gph", replace
graph export "Graphs\Calibration_men_combined.tif", name(calib_men) replace		
	
	
	
graph combine "Graphs\calibplot_newfh_prs_women.gph" ///
	"Graphs\calibplot_new_women.gph" "Graphs\calibplot_newfh_prs_men.gph" ///
	"Graphs\calibplot_new_men.gph" , cols(2) ysize(6) xsize(6) iscale(0.6) ///
	graphregion(margin(small)) name(calib_jama, replace)

graph save calib_jama "Graphs\Calibration_combined_jama.gph", replace
graph export "Graphs\Calibration_combined_jama.tif", name(calib_jama) replace		
	
	

*** Nelson-Aalen cumulative hazards ***

bysort sex: tab crc_av_10yr_quint age_gp5
bysort sex: tab crc_fh_10yr_quint age_gp5
bysort sex: tab crc_fh_prs_10yr_quint age_gp5
bysort sex: tab crc_newfh_prs_10yr_quint age_gp5
bysort sex: tab crc_new_10yr_quint age_gp5

* women
sts list if sex==0, by(crc_av_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==0, by(crc_fh_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==0, by(crc_fh_prs_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==0, by(crc_newfh_prs_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==0, by(crc_new_10yr_quint) risktable(55 65 75) cumhaz

* men
sts list if sex==1, by(crc_av_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==1, by(crc_fh_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==1, by(crc_fh_prs_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==1, by(crc_newfh_prs_10yr_quint) risktable(55 65 75) cumhaz
sts list if sex==1, by(crc_new_10yr_quint) risktable(55 65 75) cumhaz

* women - average
colorpalette tab Hue circle
sts graph if sex==0 & n_eid!=5395951, by(crc_av_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("Average risks: women", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4)	name(NA_av_10yr_f, replace) ///
	text(0.0775 32 "{bf:A}", size(*1.55) just(left))
graph save "NA_av_10yr_f" "Graphs\NA_av_10yr_women.gph", replace
graph export "Graphs\NA_av_10yr_women.tif", name(NA_av_10yr_f) replace	

		
* women - FH model 
colorpalette tab Hue circle
sts graph if sex==0 & n_eid!=5395951, by(crc_fh_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("Family history alone model: women", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_fh_10yr_f, replace) ///
	text(0.0775 32 "{bf:B}", size(*1.55) just(left))
graph save "NA_fh_10yr_f" "Graphs\NA_fh_10yr_women.gph", replace
graph export "Graphs\NA_fh_10yr_women.tif", name(NA_fh_10yr_f) replace	
	
* women - FH and PRS model 
colorpalette tab Hue circle
sts graph if sex==0 & n_eid!=5395951, by(crc_fh_prs_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("Current family history and PRS model: women", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_fh_prs_10yr_f, replace) ///
	text(0.0775 32 "{bf:C}", size(*1.55) just(left))
graph save "NA_fh_prs_10yr_f" "Graphs\NA_fh_prs_10yr_women.gph", replace
graph export "Graphs\NA_fh_prs_10yr_women.tif", name(NA_fh_prs_10yr_f) replace	

* women - new FH and PRS model
colorpalette tab Hue circle
sts graph if sex==0 & n_eid!=5395951, by(crc_newfh_prs_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("New family history and PRS model: women", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_newfh_prs_10yr_f, replace) ///
	text(0.0775 32 "{bf:D}", size(*1.55) just(left))
graph save "NA_newfh_prs_10yr_f" "Graphs\NA_newfh_prs_10yr_women.gph", replace
graph export "Graphs\NA_newfh_prs_10yr_women.tif", name(NA_newfh_prs_10yr_f) replace	

* women - new model
colorpalette tab Hue circle
sts graph if sex==0 & n_eid!=5395951 & n_eid!=4664840, ///
	by(crc_new_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("New multivariable model: women", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_new_10yr_f, replace) ///
	text(0.0775 32 "{bf:E}", size(*1.55) just(left))
graph save "NA_new_10yr_f" "Graphs\NA_new_10yr_women.gph", replace
graph export "Graphs\NA_new_10yr_women.tif", name(NA_new_10yr_f) replace	
	
graph combine "Graphs\NA_av_10yr_women.gph" "Graphs\NA_fh_10yr_women.gph" ///
	"Graphs\NA_fh_prs_10yr_women.gph" "Graphs\NA_newfh_prs_10yr_women.gph" ///
	"Graphs\NA_new_10yr_women.gph", cols(2) iscale(0.4) ///
	xsize(11) ysize(12) graphregion(margin(small)) name(NA_combined_women, replace)

graph save NA_combined_women "Graphs\NA_combined_women.gph", replace
graph export "Graphs\NA_combined_women.tif", name(NA_combined_women) replace
	
* men - average
colorpalette tab Hue circle
sts graph if sex==1 & n_eid!=1191601 & n_eid!=2719058, by(crc_av_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("Average risks: men", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_av_10yr_m, replace) ///
	text(0.0775 32 "{bf:A}", size(*1.55) just(left))
graph save "NA_av_10yr_m" "Graphs\NA_av_10yr_men.gph", replace
graph export "Graphs\NA_av_10yr_men.tif", name(NA_av_10yr_m) replace	

	
* men - FH model
colorpalette tab Hue circle
sts graph if sex==1 & n_eid!=1191601, by(crc_fh_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("Family history alone model: men", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_fh_10yr_m, replace) ///
	text(0.0775 32 "{bf:B}", size(*1.55) just(left))
graph save "NA_fh_10yr_m" "Graphs\NA_fh_10yr_men.gph", replace
graph export "Graphs\NA_fh_10yr_men.tif", name(NA_fh_10yr_m) replace	
	
* men - FH and PRS model
colorpalette tab Hue circle
sts graph if sex==1 & n_eid!=1191601, by(crc_fh_prs_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("Current family history and PRS model: men", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_fh_prs_10yr_m, replace) ///
	text(0.0775 32 "{bf:C}", size(*1.55) just(left))
graph save "NA_fh_prs_10yr_m" "Graphs\NA_fh_prs_10yr_men.gph", replace
graph export "Graphs\NA_fh_prs_10yr_men.tif", name(NA_fh_prs_10yr_m) replace	

* men - new FH and PRS model
colorpalette tab Hue circle
sts graph if sex==1 , by(crc_newfh_prs_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("New family history and PRS model: men", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_newfh_prs_10yr_m, replace) ///
	text(0.0775 32 "{bf:D}", size(*1.55) just(left))
graph save "NA_newfh_prs_10yr_m" "Graphs\NA_newfh_prs_10yr_men.gph", replace
graph export "Graphs\NA_newfh_prs_10yr_men.tif", name(NA_newfh_prs_10yr_m) replace	
	
* men - new model 
colorpalette tab Hue circle
sts graph if sex==1 & n_eid!=1135476, by(crc_new_10yr_quint) cumhaz tmin(40) ///
	plot1opts(lcolor("`r(p19)'") lwidth(thick)) /// use medthick for single graph
	plot2opts(lcolor("`r(p3)'") lwidth(thick)) ///
	plot3opts(lcolor("`r(p6)'") lwidth(thick)) ///
	plot4opts(lcolor("`r(p9)'") lwidth(thick)) ///
	plot5opts(lcolor("`r(p13)'") lwidth(thick)) ///
	yscale(range(0 0.07))  ///
	ylabel(0(0.01)0.07, format(%5.2f) nogrid angle(horizontal)) ///
	xscale(range(39 81)) xlabel(40(5)80) ///
	title("New multivariable model: men", size(*1)) ///
	ytitle("Nelson–Aalen cumulative hazard", height(4)) ///
	xtitle("Age (years)", height(4)) ///
	legend(label(1 "Quintile 1") label(2 "Quintile 2") ///
		label(3 "Quintile 3") label(4 "Quintile 4") label(5 "Quintile 5") ///
		order(5 4 3 2 1) size(*0.9) symxsize(7) keygap(1.5) rowgap(0) ///
		cols(1) ring(0) bplacement(nwest) region(lwidth(none))) ///
	xsize(5.5) ysize(4) name(NA_new_10yr_m, replace) ///
	text(0.0775 32 "{bf:E}", size(*1.55) just(left))
graph save "NA_new_10yr_m" "Graphs\NA_new_10yr_men.gph", replace
graph export "Graphs\NA_new_10yr_men.tif", name(NA_new_10yr_m) replace	


graph combine "Graphs\NA_av_10yr_men.gph" "Graphs\NA_fh_10yr_men.gph" ///
	"Graphs\NA_fh_prs_10yr_men.gph" "Graphs\NA_newfh_prs_10yr_men.gph" ///
	"Graphs\NA_new_10yr_men.gph", cols(2) iscale(0.4) ///
	xsize(11) ysize(12) graphregion(margin(small)) name(NA_combined_men, replace)

graph save NA_combined_men "Graphs\NA_combined_men.gph", replace
graph export "Graphs\NA_combined_men.tif", name(NA_combined_men) replace

* 10-year risk distributions

colorpalette tab Hue circle	
twoway histogram crc_av_10yr if sex==0, start(0) width(0.002) recast(line) ///
		lcolor("`r(p19)'") lwidth(medthick) || ///
	histogram crc_fh_10yr if sex==0, start(0) width(0.002) recast(line) ///
		lcolor("`r(p3)'") lwidth(medthick) || ///
	histogram crc_fh_prs_10yr if sex==0 & crc_fh_prs_10yr<0.06, start(0) ///
		width(0.002) recast(line) lcolor("`r(p6)'") lwidth(medthick) || ///
	histogram crc_newfh_prs_10yr if sex==0 & crc_newfh_prs_10yr<0.06, start(0) ///
		width(0.002) recast(line) lcolor("`r(p9)'") lwidth(medthick) || ///
	histogram crc_new_10yr if sex==0 & crc_new_10yr<0.06, start(0) width(0.002) ///
		recast(line) lcolor("`r(p13)'") lwidth(medthick) ///
	xscale(range(0 0.06)) xlabel(0(0.01)0.06, format(%4.2f)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5" 0.06 "6") ///
	yscale(range(0 80)) ylabel(0(10)80, angle(horizontal) format(%3.0f)) ///
	title("Colorectal cancer: women") xtitle("10-year risk (%)", height(5)) ///
	ytitle(, height(5)) ///
	legend(region(lwidth(none)) order(1 2 3 4 5) cols(2) size(*0.75) ///
	label(1 "Average risks") label(2 "Family history alone") ///
	label(3 "Current family history and PRS") ///
	label(4 "New family history and PRS") label(5 "New multivariable")) ///
	name(risk_dist_f, replace)

graph save "risk_dist_f" "Graphs\Risk_dist_10yr_women.gph", replace
graph export "Graphs\Risk_dist_10yr_women.tif", name(risk_dist_f) replace		
	
colorpalette tab Hue circle	
twoway histogram crc_av_10yr if sex==1, start(0) width(0.002) recast(line) ///
		lcolor("`r(p19)'") lwidth(medthick) || ///
	histogram crc_fh_10yr if sex==1, start(0) width(0.002) recast(line) ///
		lcolor("`r(p3)'") lwidth(medthick) || ///
	histogram crc_fh_prs_10yr if sex==1 & crc_fh_prs_10yr<0.06, start(0) ///
		width(0.002) recast(line) lcolor("`r(p6)'") lwidth(medthick) || ///
	histogram crc_newfh_prs_10yr if sex==1 & crc_newfh_prs_10yr<0.06, start(0) ///
		width(0.002) recast(line) lcolor("`r(p9)'") lwidth(medthick) || ///
	histogram crc_new_10yr if sex==1 & crc_new_10yr<0.06, start(0) width(0.002) ///
		recast(line) lcolor("`r(p13)'") lwidth(medthick) ///
	xscale(range(0 0.06)) xlabel(0(0.01)0.06, format(%4.2f)) ///
	xlabel(0 "0" 0.01 "1" 0.02 "2" 0.03 "3" 0.04 "4" 0.05 "5" 0.06 "6") ///
	yscale(range(0 80)) ylabel(0(10)80, angle(horizontal) format(%3.0f)) ///
	title("Colorectal cancer: men") xtitle("10-year risk (%)", height(5)) ///
	ytitle(, height(5)) ///
	legend(region(lwidth(none)) order(1 2 3 4 5) cols(2) size(*0.75) ///
	label(1 "Average risks") label(2 "Family history alone") ///
	label(3 "Current family history and PRS") ///
	label(4 "New family history and PRS") label(5 "New multivariable")) ///
	name(risk_dist_m, replace)

graph save "risk_dist_m" "Graphs\Risk_dist_10yr_men.gph", replace
graph export "Graphs\Risk_dist_10yr_men.tif", name(risk_dist_m) replace		
	
* exclude outlier with crc_fh_prs_10yr>0.2 (n_eid=1579709)	

*** scatterplots

* women
colorpalette tab Hue circle
twoway scatter crc_fh_prs_10yr crc_new_10yr if cr10==0 & sex==1,  ///
	msize(tiny) mcolor("`r(p19)'") msymbol(o) || ///
	scatter crc_fh_prs_10yr crc_new_10yr if cr10==1 & sex==1 & ///
		crc_fh_prs_10yr<0.2, aspect(1) ///
	msize(tiny) mcolor("`r(p13)'") msymbol(o) ///
	xscale(range(0 0.2)) xlabel(0 "0" 0.05 "5" 0.1 "10" 0.15 "15" 0.2 "20") ///
	yscale(range(0 0.2)) ylabel(0 "0" 0.05 "5" 0.1 "10" 0.15 "15" 0.2 "20", ///
		angle(horizontal)) ///
	xsize(4) ysize(4.15) ///
	xtitle("New multivariable model: 10-year risk (%)", height(5)) ///
	ytitle("Current family history and PRS model:" "10-year risk (%)", height(7)) ///
	title("Women") ///
	legend(region(lwidth(none)) size(*0.75) col(1) ring(0) position(5) ///
		label(1 "Unaffected") label(2 "Affected")) ///
	name(scatter_new_newfhprs_f, replace)
	
graph save "scatter_new_newfhprs_f" "Graphs\Scatter_new_newfhprs_f.gph", replace
graph export "Graphs\Scatter_new_newfhprs_f.tif", name(scatter_new_newfhprs_f) replace

* men
colorpalette tab Hue circle
twoway scatter crc_fh_prs_10yr crc_new_10yr if cr10==0 & sex==1,  ///
	msize(tiny) mcolor("`r(p19)'") msymbol(o) || ///
	scatter crc_fh_prs_10yr crc_new_10yr if cr10==1 & sex==1 & ///
		crc_fh_prs_10yr<0.2, aspect(1) ///
	msize(tiny) mcolor("`r(p13)'") msymbol(o) ///
	xscale(range(0 0.2)) xlabel(0 "0" 0.05 "5" 0.1 "10" 0.15 "15" 0.2 "20") ///
	yscale(range(0 0.2)) ylabel(0 "0" 0.05 "5" 0.1 "10" 0.15 "15" 0.2 "20", ///
		angle(horizontal)) ///
	xsize(4) ysize(4.15) ///
	xtitle("New multivariable model: 10-year risk (%)", height(5)) ///
	ytitle("Current family history and PRS model:" "10-year risk (%)", height(7)) ///
	title("Men") ///
	legend(region(lwidth(none)) size(*0.75) col(1) ring(0) position(5) ///
		label(1 "Unaffected") label(2 "Affected")) ///
	name(scatter_new_newfhprs_m, replace)

graph save "scatter_new_newfhprs_m" "Graphs\Scatter_new_newfhprs_m.gph", replace
graph export "Graphs\Scatter_new_newfhprs_m.tif", name(scatter_new_newfhprs_m) replace
	
	
*** histograms

* all - FH alone
colorpalette tab Hue circle
twoway histogram crc_fh_10yr if cr10==0 & sex==0, start(0) width(0.002) ///
			recast(line) lcolor("`r(p13)'") lwidth(medthick) || ///
		histogram crc_fh_10yr if cr10==1 & sex==0, start(0) width(0.005) ///
			recast(line) lcolor("`r(p13)'") lpattern(dash) lwidth(medthick) || ///
	histogram crc_fh_10yr if cr10==0 & sex==1, start(0) width(0.002) ///
		recast(line) lcolor("`r(p19)'") lwidth(medthick) || ///
	histogram crc_fh_10yr if cr10==1 & sex==1, start(0) ///
		width(0.005) recast(line) lcolor("`r(p19)'") lpattern(dash) lwidth(medthick) ///
	xscale(range(0 0.15)) xlabel(0 "0" 0.02 "2" 0.04 "4" 0.06 "6" 0.08 "8" ///
		0.1 "10" 0.12 "12" 0.14 "14") ///
	ylabel(0(10)80, angle(horizontal)) ///
	title("Family history alone model") ///
	xtitle("10-year risk (%)", height(5)) ///
	legend(region(lwidth(none)) size(*0.75) rowgap(*0.2) symxsize(5) ///
		col(1) ring(0) position(1) ///
		label(1 "Unaffected women") label(2 "Affected women") ///
		label(3 "Unaffected men") label(4 "Affected men")) ///
	name(histogram_fh, replace)

graph save "histogram_fh" "Graphs\histogram_fh.gph", replace
graph export "Graphs\Histogram_fh.tif", name(histogram_fh) replace	
	
* all - current FH and PRS
colorpalette tab Hue circle
twoway histogram crc_fh_prs_10yr if cr10==0 & sex==0, start(0) width(0.002) ///
			recast(line) lcolor("`r(p13)'") lwidth(medthick) || ///
		histogram crc_fh_prs_10yr if cr10==1 & sex==0, start(0) width(0.005) ///
			recast(line) lcolor("`r(p13)'") lpattern(dash) lwidth(medthick) || ///
	histogram crc_fh_prs_10yr if cr10==0 & sex==1, start(0) width(0.002) ///
		recast(line) lcolor("`r(p19)'") lwidth(medthick) || ///
	histogram crc_fh_prs_10yr if cr10==1 & sex==1, start(0) ///
		width(0.005) recast(line) lcolor("`r(p19)'") lpattern(dash) lwidth(medthick) ///
	xscale(range(0 0.15)) xlabel(0 "0" 0.02 "2" 0.04 "4" 0.06 "6" 0.08 "8" ///
		0.1 "10" 0.12 "12" 0.14 "14") ///
	ylabel(0(10)80, angle(horizontal)) ///
	title("Current family history and PRS model") ///
	xtitle("10-year risk (%)", height(5)) ///
	legend(region(lwidth(none)) size(*0.75) rowgap(*0.2) symxsize(5) ///
		col(1) ring(0) position(1) ///
		label(1 "Unaffected women") label(2 "Affected women") ///
		label(3 "Unaffected men") label(4 "Affected men")) ///
	name(histogram_fh_prs, replace)

graph save "histogram_fh_prs" "Graphs\histogram_fh_prs.gph", replace
graph export "Graphs\Histogram_fh_prs.tif", name(histogram_fh_prs) replace	
	
	
* all - new FH and PRS
colorpalette tab Hue circle
twoway histogram crc_newfh_prs_10yr if cr10==0 & sex==0, start(0) width(0.002) ///
			recast(line) lcolor("`r(p13)'") lwidth(medthick) || ///
		histogram crc_newfh_prs_10yr if cr10==1 & sex==0, start(0) width(0.005) ///
			recast(line) lcolor("`r(p13)'") lpattern(dash) lwidth(medthick) || ///
	histogram crc_newfh_prs_10yr if cr10==0 & sex==1, start(0) width(0.002) ///
		recast(line) lcolor("`r(p19)'") lwidth(medthick) || ///
	histogram crc_newfh_prs_10yr if cr10==1 & sex==1, start(0) ///
		width(0.005) recast(line) lcolor("`r(p19)'") lpattern(dash) lwidth(medthick) ///
	xscale(range(0 0.15)) xlabel(0 "0" 0.02 "2" 0.04 "4" 0.06 "6" 0.08 "8" ///
		0.1 "10" 0.12 "12" 0.14 "14") ///
	ylabel(0(10)80, angle(horizontal)) ///
	title("New family history and PRS model") ///
	xtitle("10-year risk (%)", height(5)) ///
	legend(region(lwidth(none)) size(*0.75) rowgap(*0.2) symxsize(5) ///
		col(1) ring(0) position(1) ///
		label(1 "Unaffected women") label(2 "Affected women") ///
		label(3 "Unaffected men") label(4 "Affected men")) ///
	name(histogram_newfh_prs, replace)
	
graph save "histogram_newfh_prs" "Graphs\histogram_newfh_prs.gph", replace
graph export "Graphs\Histogram_newfh_prs.tif", name(histogram_newfh_prs) replace	
	
* all - new multivariable
colorpalette tab Hue circle
twoway histogram crc_new_10yr if cr10==0 & sex==0, start(0) width(0.002) ///
			recast(line) lcolor("`r(p13)'") lwidth(medthick) || ///
		histogram crc_new_10yr if cr10==1 & sex==0, start(0) width(0.005) ///
			recast(line) lcolor("`r(p13)'") lpattern(dash) lwidth(medthick) || ///
	histogram crc_new_10yr if cr10==0 & sex==1, start(0) width(0.002) ///
		recast(line) lcolor("`r(p19)'") lwidth(medthick) || ///
	histogram crc_new_10yr if cr10==1 & sex==1, start(0) ///
		width(0.005) recast(line) lcolor("`r(p19)'") lpattern(dash) lwidth(medthick) ///
	xscale(range(0 0.15)) xlabel(0 "0" 0.02 "2" 0.04 "4" 0.06 "6" 0.08 "8" ///
		0.1 "10" 0.12 "12" 0.14 "14") ///
	ylabel(0(10)80, angle(horizontal)) ///
	title("New multivariable model") ///
	xtitle("10-year risk (%)", height(5)) ///
	legend(region(lwidth(none)) size(*0.75) rowgap(*0.2) symxsize(5) ///
		col(1) ring(0) position(1) ///
		label(1 "Unaffected women") label(2 "Affected women") ///
		label(3 "Unaffected men") label(4 "Affected men")) ///
	name(histogram_new, replace)

graph save "histogram_new" "Graphs\histogram_new.gph", replace
graph export "Graphs\Histogram_new.tif", name(histogram_new) replace	

	
	
log close testing
exit