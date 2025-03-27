*** UK Biobank -- colorectal cancer modelling PRS using imputed data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Model_FH_imputed_`filedate'.log", ///
	name(model_prs_impute) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m


*** PRS MODEL - Women ***

mi estimate, hr saving(models\model_prs_f, replace): stcox std_prs_xb if sex==0
mi estimate, nohr

* check for TVC
mi estimate, hr: stcox std_prs_xb if sex==0, tvc(std_prs_xb)

stcox std_prs_xb if sex==0 & _mi_m==1, hr
estat phtest, detail


*** PRS  MODEL - Men ***

mi estimate, hr saving(models\model_prs_m, replace): stcox std_prs_xb if sex==1
mi estimate, nohr

* check for TVC
mi estimate, hr: stcox std_prs_xb if sex==1, tvc(std_prs_xb)

stcox std_prs_xb if sex==1 & _mi_m==1, hr
estat phtest, detail

estat phtest, plot(std_prs_xb) name(std_prs_xb, replace) // don't need interaction


*** model fit ***

mi extract 1

* women
preserve
keep if sex==0

stcox std_prs_xb, nohr mgale(mg)

predict cs, csnell // Cox-Snell residuals

drop mg

stset cs, failure(cr==1)
sts generate H = na // Nelson-Aalen cumulative hazard
twoway line H cs if cs<=0.05, sort lcolor(orange) || ///
	line cs cs if cs<=.05, sort lcolor(ebblue) ///
	aspect(1) xsize(3.5) xtitle(, margin(t +2)) ///
	xscale(range(0 0.05)) xlabels(0(0.01)0.05, format(%03.2f)) ///
	yscale(range(0 0.05)) ylabels(0(0.01)0.05, format(%03.2f) angle(h))  ///
	legend(cols(1) size(small) region(lwidth(none))) ///
		name(NA_newprs_f, replace)
graph save "NA_newprs_f" "Graphs\CS_NA_newprs_f.gph", replace
graph export "Graphs\CS_NA_newprs_f.tif", as(tif) ///
	name("NA_newprs_f") replace

restore

* men
preserve
keep if sex==1

stcox std_prs_xb, nohr mgale(mg)

predict cs, csnell // Cox-Snell residuals

drop mg

stset cs, failure(cr==1)
sts generate H = na // Nelson-Aalen cumulative hazard
twoway line H cs if cs<=0.15, sort lcolor(orange) || ///
	line cs cs if cs<=.15, sort lcolor(ebblue) ///
	aspect(1) xsize(3.5) xtitle(, margin(t +2)) ///
	xscale(range(0 0.14)) xlabels(0(0.02)0.14, format(%03.2f)) ///
	yscale(range(0 0.14)) ylabels(0(0.02)0.14, format(%03.2f) angle(h))  ///
	legend(cols(1) size(small) region(lwidth(none))) ///
		name(NA_newprs_m, replace)
graph save "NA_newprs_m" "Graphs\CS_NA_newprs_m.gph", replace
graph export "Graphs\CS_NA_newprs_m.tif", as(tif) ///
	name("NA_newprs_m") replace

restore

log close model_prs_impute
