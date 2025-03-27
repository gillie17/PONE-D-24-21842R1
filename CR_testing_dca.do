*** DECISION CURVE ANALYSIS ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Testing_DCA_`filedate'.log", ///
	name(testing) replace

use "CR_testing_analysis.dta", clear

keep if ancestry_code==9

* use time-to-event data for this analysis (fup10_sir_time)

* competing risk variable
generate status = 0
replace status = 1 if cr10==1
replace status = 2 if cr10==0 & fup10_sir_time<10

stset fup10_sir_time, id(n_eid) failure(status==1) // limit to 10 years of follow-up

* competing risk adjustment -- compet1(2) -- doesn't work for some reason
	
* women
colorpalette tab Hue circle
stdca crc_fh_10yr crc_newfh_prs_10yr crc_new_10yr ///
	if sex==0, timepoint(10) smooth xstart(0.001) xstop(0.06) ymin(-0.008) ///
	lcolor("`r(p19)'" "`r(p3)'" "`r(p6)'" "`r(p9)'" "`r(p13)'") ///
	lwidth(medthick medthick medthick medthick medthick medthick) ///
	xline(0.01, lwidth(vthin) lcolor(gs8)) ///
	xline(0.02, lwidth(vthin) lcolor(gs8)) ///
	xlabel(0 "0.0" 0.01 "1.0" 0.02 "2.0" 0.03 "3.0" 0.04 "4.0" 0.05 "5.0") ///
	yscale(range(-0.001 0.008)) ///
	ylabel(-0.002(0.002)0.008, angle(h) format("%9.3f")) ///
	xtitle("Threshold probability (%)") ///
	ytitle("Net benefit", height(0)) ///
	title("Decision curves: women", size(0.75*)) ///
	legend(size(small) cols(1) region(lwidth(none)) ///
	symxsize(7) keygap(1.5) rowgap(*0.25) label(1 "Screen all") ///
		label(2 "Screen none")  label(3 "Family history alone") ///
		label(4 "New family history and PRS") ///
		label(5 "New multivariable") ring(0) bplacement(neast)) ///
	text(0.009 -0.01 "{bf:A}", size(*1.5) just(left)) ///
	name(DCA_f, replace)

graph save "DCA_f" "Graphs\DCA_women.gph", replace
graph export "Graphs\DCA_women.tif", as(tif) name("DCA_f") replace


* men	
colorpalette tab Hue circle
stdca crc_fh_10yr crc_newfh_prs_10yr crc_new_10yr ///
	if sex==1, timepoint(10) smooth xstart(0.001) xstop(0.06) ymin(-0.01) ///
	lcolor("`r(p19)'" "`r(p3)'" "`r(p6)'" "`r(p9)'" "`r(p13)'") ///
	lwidth(medthick medthick medthick medthick medthick medthick) ///
	xline(0.01, lwidth(vthin) lcolor(gs8)) ///
	xline(0.02, lwidth(vthin) lcolor(gs8)) ///
	xlabel(0 "0.0" 0.01 "1.0" 0.02 "2.0" 0.03 "3.0" 0.04 "4.0" 0.05 "5.0") ///
	yscale(range(-0.005 0.015)) ///
	ylabel(-0.005(0.005)0.015, angle(h) format("%9.3f")) ///
	xtitle("Threshold probability (%)") ///
	ytitle("Net benefit", height(0)) ///	
	title("Decision curves: men", size(0.75*)) ///
	legend(size(small) cols(1) region(lwidth(none)) ///
	symxsize(7) keygap(1.5) rowgap(*0.25) label(1 "Screen all") ///
		label(2 "Screen none")  label(3 "Family history alone") ///
		label(4 "New family history and PRS") ///
		label(5 "New multivariable") ring(0) bplacement(neast)) ///
	text(0.017 -0.01 "{bf:B}", size(*1.5) just(left)) ///
	name(DCA_m, replace)

graph save "DCA_m" "Graphs\DCA_men.gph", replace
graph export "Graphs\DCA_men.tif", as(tif) name("DCA_m") replace


graph combine "Graphs\DCA_women.gph" "Graphs\DCA_men.gph", cols(1) iscale(0.8) ///
	xsize(5.5) ysize(8) graphregion(margin(none)) name(DCA_combined, replace)

graph save DCA_combined "Graphs\DCA_combined.gph", replace
graph export "Graphs\DCA_combined.tif", name(DCA_combined) replace

log close testing
exit
