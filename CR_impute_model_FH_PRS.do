*** UK Biobank -- colorectal cancer modelling FH and PRS using imputed data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Model_FH_PRS_imputed_`filedate'.log", ///
	name(model_fh_prs_impute) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m


*** FH AND PRS MODEL - Women ***

mi estimate, hr saving(models\model_fh_prs_f, replace): stcox std_prs_xb ///
	cr_deg1 if sex==0
mi estimate, nohr

* check for TVC
mi estimate, hr: stcox std_prs_xb cr_deg1 if sex==0, tvc(std_prs_xb cr_deg1)

stcox std_prs_xb cr_deg1 if sex==0 & _mi_m==1, hr
estat phtest, detail


*** FH AND PRS MODEL - Men ***

mi estimate, hr saving(models\model_fh_prs_m, replace): stcox std_prs_xb ///
	cr_deg1 if sex==1
mi estimate, nohr

* check for TVC
mi estimate, hr: stcox std_prs_xb cr_deg1 if sex==1, tvc(std_prs_xb cr_deg1)

stcox std_prs_xb cr_deg1 if sex==1 & _mi_m==1, hr
estat phtest, detail

colorpalette tab Hue circle
estat phtest, plot(std_prs_xb)  /// don't need interaction
	mcolor("`r(p13)'") msize(vsmall) lineopts(lcolor("`r(p19)'") lwidth(thick)) ///
	xlabel(40(5)80) ///
	ylabel(-4(1)4, angle(h)) ///
	ytitle("Schoenfeld residuals", height(4)) note("") ///
	xtitle("Age (years)", height(5)) ///
	title("New family history and PRS model for men:" "140-SNP PRS", ///
		size(0.75*)) ///
	text(5.25 36 "{bf:A}", size(*1.5) just(left)) ///
	name(std_prs_xb_fhprs, replace)

graph save "std_prs_xb_fhprs" "Graphs\Schoenfeld_prs_m_fhprs.gph", replace
graph export "Graphs\Schoenfeld_prs_m_fhprs.tif", as(tif) ///
	name("std_prs_xb_fhprs") replace

colorpalette tab Hue circle
estat phtest, plot(cr_deg1) /// don't need interaction
	mcolor("`r(p13)'") msize(vsmall) lineopts(lcolor("`r(p19)'") lwidth(thick)) ///
	xlabel(40(5)80) ///
	ylabel(-2(2)8, angle(h)) ///
	ytitle("Schoenfeld residuals", height(4)) note("") ///
	xtitle("Age (years)", height(5)) ///
	title("New family history and PRS model for men:" "first-degree family history", ///
		size(0.75*)) ///
	text(9.45 36 "{bf:B}", size(*1.5) just(left)) ///
	name(cr_deg1_fhprs, replace)

graph save "cr_deg1_fhprs" "Graphs\Schoenfeld_deg1_m_fhprs.gph", replace
graph export "Graphs\Schoenfeld_deg1_m_fhprs.tif", as(tif) ///
	name("cr_deg1_fhprs") replace



*** model fit ***

mi extract 1

* women
preserve
keep if sex==0

stcox std_prs_xb cr_deg1, nohr mgale(mg)

predict cs, csnell // Cox-Snell residuals

drop mg

stset cs, failure(cr==1)
sts generate H = na // Nelson-Aalen cumulative hazard

colorpalette tab Hue circle
twoway line cs cs if cs<=.05, sort lcolor("`r(p19)'") lwidth(thick) || ///
	line H cs if cs<=0.05, sort lcolor("`r(p13)'") lwidth(thick) ///
	aspect(1) xsize(3.75) ysize(4) xtitle(, margin(t +2)) ///
	xscale(range(0 0.05)) xlabels(0(0.01)0.05, format(%03.2f)) ///
	yscale(range(0 0.05)) ylabels(0(0.01)0.05, format(%03.2f) angle(h)) ///
	title("New family history and PRS model: women", size(*0.75)) ///
	legend(cols(1) size(small) region(lwidth(none)) ring(0) position(5) ///
	symxsize(7) rowgap(*0.25)) ///
	text(0.055 -0.007 "{bf:A}", size(*1.5) just(left)) ///
	name(NA_newfh_prs_f, replace)
graph save "NA_newfh_prs_f" "Graphs\CS_NA_newfh_prs_f.gph", replace
graph export "Graphs\CS_NA_newfh_prs_f.tif", as(tif) ///
	name("NA_newfh_prs_f") replace

restore

* men
preserve
keep if sex==1

stcox std_prs_xb cr_deg1, nohr mgale(mg)

predict cs, csnell // Cox-Snell residuals

drop mg

stset cs, failure(cr==1)
sts generate H = na // Nelson-Aalen cumulative hazard

colorpalette tab Hue circle
twoway line cs cs if cs<=.15, sort lcolor("`r(p19)'") lwidth(thick) || ///
	line H cs if cs<=0.15, sort lcolor("`r(p13)'") lwidth(thick) ///
	aspect(1) xsize(3.75) ysize(4) xtitle(, margin(t +2)) ///
	xscale(range(0 0.14)) xlabels(0(0.02)0.14, format(%03.2f)) ///
	yscale(range(0 0.14)) ylabels(0(0.02)0.14, format(%03.2f) angle(h)) ///
	title("New family history and PRS model: men", size(*0.75)) ///
	legend(cols(1) size(small) region(lwidth(none)) ring(0) position(5) ///
		symxsize(7) rowgap(*0.25)) ///
	text(0.154 -0.019 "{bf:B}", size(*1.5) just(left)) ///
	name(NA_newfh_prs_m, replace)
graph save "NA_newfh_prs_m" "Graphs\CS_NA_newfh_prs_m.gph", replace
graph export "Graphs\CS_NA_newfh_prs_m.tif", as(tif) ///
	name("NA_newfh_prs_m") replace

restore

log close model_fh_prs_impute


