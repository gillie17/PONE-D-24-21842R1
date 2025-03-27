*** UK Biobank -- colorectal cancer Kevin's models using imputed data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Model_Kevin_imputed_`filedate'.log", ///
	name(model_kevin_impute) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m


*** Kevin's model - women ***

mi estimate, saving(models\model_kevin_f, replace): stcox std_prs_xb ///
	cr_deg1 smoke_ever screen10_new n_30870_0_0_c if sex==0
mi estimate, hr base
mi estimate, coeflegend 

how_many_imputations, cv_se(0.05) 

* check for TVC
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever screen10_new ///
	n_30870_0_0_c if sex==0, tvc(std_prs_xb cr_deg1 smoke_ever screen10_new ///
	n_30870_0_0_c)	

* proportional hazards test in first imputation dataset - no problems

stcox std_prs_xb cr_deg1 smoke_ever screen10_new n_30870_0_0_c ///
	if sex==0 & _mi_m==1, hr
estat phtest, detail


*** model fit ***

mi extract 1

preserve
keep if sex==0

stcox std_prs_xb cr_deg1 smoke_ever screen10_new n_30870_0_0_c, nohr mgale(mg)

predict cs, csnell // Cox-Snell residuals
drop mg

stset cs, failure(cr==1)
sts generate H = na // Nelson-Aalen cumulative hazard

colorpalette tab Hue circle
twoway line cs cs if cs<=.05, sort lcolor("`r(p19)'") lwidth(thick) || ///
	line H cs if cs<=0.05, sort lcolor("`r(p13)'") lwidth(thick)  ///
		aspect(1) xsize(3.75) ysize(4) xtitle(, margin(t +2)) ///
	xscale(range(0 0.05)) xlabels(0(0.01)0.05, format(%03.2f)) ///
	yscale(range(0 0.05)) ylabels(0(0.01)0.05, format(%03.2f) angle(h)) ///
	title("New multivariable model: women", size(*0.75)) ///
	legend(cols(1) size(*0.8) region(lwidth(none)) ring(0) position(5) ///
	symxsize(7) rowgap(*0.25)) ///
	text(0.055 -0.007 "{bf:C}", size(*1.5) just(left)) ///
	name(NA_kevin_f, replace)
graph save "NA_kevin_f" "Graphs\CS_NA_kevin_f.gph", ///
	replace
graph export "Graphs\CS_NA_kevin_f.tif", as(tif) ///
	name("NA_kevin_f") replace

restore


*** Kevin's model - men ***

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m


mi estimate, saving(models\model_kevin_m, replace): stcox std_prs_xb ///
	cr_deg1 smoke_ever screen10_new ln_bmi_c if sex==1
mi estimate, hr base
mi estimate, coeflegend 

how_many_imputations, cv_se(0.05) 

* check for TVC
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever screen10_new ln_bmi_c ///
	if sex==1, tvc(std_prs_xb cr_deg1 smoke_ever screen10_new ln_bmi_c )

stcox std_prs_xb cr_deg1 smoke_ever screen10_new ln_bmi_c if sex==1 & ///
	_mi_m==1, hr
estat phtest, detail

colorpalette tab Hue circle
estat phtest, plot(std_prs_xb)  /// don't need interaction
	mcolor("`r(p13)'") msize(vsmall) lineopts(lcolor("`r(p19)'") lwidth(thick)) ///
	xlabel(40(5)80) ///
	ylabel(-4(1)4, angle(h)) ///
	ytitle("Schoenfeld residuals", height(4)) note("") ///
	xtitle("Age (years)", height(5)) ///
	title("New multivariable model for men:" "140-SNP PRS", ///
		size(0.75*)) ///
	text(5.25 36 "{bf:C}", size(*1.5) just(left)) ///
	name(std_prs_xb_new, replace)

graph save "std_prs_xb_new" "Graphs\Schoenfeld_prs_m_new.gph", replace
graph export "Graphs\Schoenfeld_prs_m_new.tif", as(tif) ///
	name("std_prs_xb_new") replace

colorpalette tab Hue circle
estat phtest, plot(cr_deg1) /// don't need interaction
	mcolor("`r(p13)'") msize(vsmall) lineopts(lcolor("`r(p19)'") lwidth(thick)) ///
	xlabel(40(5)80) ///
	ylabel(-2(2)8, angle(h)) ///
	ytitle("Schoenfeld residuals", height(4)) note("") ///
	xtitle("Age (years)", height(5)) ///
	title("New multivariable model for men:" "first-degree family history", ///
		size(0.75*)) ///
	text(9.45 36 "{bf:D}", size(*1.5) just(left)) ///
	name(cr_deg1_new, replace)

graph save "cr_deg1_new" "Graphs\Schoenfeld_deg1_m_new.gph", replace
graph export "Graphs\Schoenfeld_deg1_m_new.tif", as(tif) ///
	name("cr_deg1_new") replace

*** model fit ***

mi extract 1

preserve
keep if sex==1

stcox std_prs_xb cr_deg1 smoke_ever screen10_new ln_bmi_c, nohr mgale(mg)

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
	title("New multivariable model: men", size(*0.75)) ///
	legend(cols(1) size(small) region(lwidth(none)) ring(0) position(5) ///
		symxsize(7) rowgap(*0.25)) ///
	text(0.1558 -0.019 "{bf:D}", size(*1.5) just(left)) ///
	name(NA_kevin_m, replace)
graph save "NA_kevin_m" "Graphs\CS_NA_kevin_m.gph", ///
	replace
graph export "Graphs\CS_NA_kevin_m.tif", as(tif) ///
	name("NA_kevin_m") replace

restore

log close model_kevin_impute


* combine graphs for paper


* 2 lines for graph titles
graph combine "Graphs\Schoenfeld_prs_m_fhprs.gph" ///
	"Graphs\Schoenfeld_deg1_m_fhprs.gph" ///
	"Graphs\Schoenfeld_prs_m_new.gph" "Graphs\Schoenfeld_deg1_m_new.gph", ///
	cols(2) iscale(0.45) xsize(8) ysize(6) graphregion(margin(small)) ///
	name(Schoenfeld_combined)

graph save "Schoenfeld_combined" "Graphs\Schoenfeld_combined.gph", replace
graph export "Graphs\Schoenfeld_combined.tif", as(tif) ///
	name("Schoenfeld_combined") replace	
	

graph combine "Graphs\CS_NA_newfh_prs_f.gph" "Graphs\CS_NA_newfh_prs_m.gph" ///
	"Graphs\CS_NA_kevin_f.gph" "Graphs\CS_NA_kevin_m.gph", ///
	cols(2) iscale(0.45) xsize(11) ysize(12) graphregion(margin(small)) ///
	name(CS_NA_combined, replace)
	
graph save "CS_NA_combined" "Graphs\CS_NA_combined.gph", replace
graph export "Graphs\CS_NA_combined.tif", as(tif) ///
	name("CS_NA_combined") replace

