*** CRC model model performance in people of UK ancestry ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\SIR_analysis_`filedate'.log", ///
	name(sir) replace

use "CR_testing_analysis.dta", clear

generate crc_av_1yr=crc_av_10yr/10
generate crc_fh_1yr=crc_fh_10yr/10
generate crc_fh_prs_1yr=crc_fh_prs_10yr/10
generate crc_newfh_prs_1yr=crc_newfh_prs_10yr/10
generate crc_new_1yr=crc_new_10yr/10

keep if ancestry_code==9


*** new groups for the graphs

generate crc_av_10yr_gp=crc_av_10yr_quint
replace crc_av_10yr_gp=6 if crc_av_10yr_dec==10

generate crc_fh_10yr_gp=crc_fh_10yr_quint
replace crc_fh_10yr_gp=6 if crc_fh_10yr_dec==10

generate crc_fh_prs_10yr_gp=crc_fh_prs_10yr_quint
replace crc_fh_prs_10yr_gp=6 if crc_fh_prs_10yr_dec==10

generate crc_newfh_prs_10yr_gp=crc_newfh_prs_10yr_quint
replace crc_newfh_prs_10yr_gp=6 if crc_newfh_prs_10yr_dec==10

generate crc_new_10yr_gp=crc_new_10yr_quint
replace crc_new_10yr_gp=6 if crc_new_10yr_dec==10

*** WOMEN ***

* average - women

tabstat crc_av_10yr if sex==0, by(crc_av_10yr_gp) ///
	statistics(p50) nototal save
matrix f_av1=r(Stat1)
matrix f_av2=r(Stat2)
matrix f_av3=r(Stat3)
matrix f_av4=r(Stat4)
matrix f_av5=r(Stat5)
matrix f_av6=r(Stat6)

scalar f_av1_med=f_av1[1,1]
scalar f_av2_med=f_av2[1,1]
scalar f_av3_med=f_av3[1,1]
scalar f_av4_med=f_av4[1,1]
scalar f_av5_med=f_av5[1,1]
scalar f_av6_med=f_av6[1,1]

replace crc_av_10yr_gp=f_av1_med if sex==0 & crc_av_10yr_gp==1
replace crc_av_10yr_gp=f_av2_med if sex==0 & crc_av_10yr_gp==2
replace crc_av_10yr_gp=f_av3_med if sex==0 & crc_av_10yr_gp==3
replace crc_av_10yr_gp=f_av4_med if sex==0 & crc_av_10yr_gp==4
replace crc_av_10yr_gp=f_av5_med if sex==0 & crc_av_10yr_gp==5
replace crc_av_10yr_gp=f_av6_med if sex==0 & crc_av_10yr_gp==6

* current FH - women

tabstat crc_fh_10yr if sex==0, by(crc_fh_10yr_gp) ///
	statistics(p50) nototal save
matrix f_fh1=r(Stat1)
matrix f_fh2=r(Stat2)
matrix f_fh3=r(Stat3)
matrix f_fh4=r(Stat4)
matrix f_fh5=r(Stat5)
matrix f_fh6=r(Stat6)

scalar f_fh1_med=f_fh1[1,1]
scalar f_fh2_med=f_fh2[1,1]
scalar f_fh3_med=f_fh3[1,1]
scalar f_fh4_med=f_fh4[1,1]
scalar f_fh5_med=f_fh5[1,1]
scalar f_fh6_med=f_fh6[1,1]

replace crc_fh_10yr_gp=f_fh1_med if sex==0 & crc_fh_10yr_gp==1
replace crc_fh_10yr_gp=f_fh2_med if sex==0 & crc_fh_10yr_gp==2
replace crc_fh_10yr_gp=f_fh3_med if sex==0 & crc_fh_10yr_gp==3
replace crc_fh_10yr_gp=f_fh4_med if sex==0 & crc_fh_10yr_gp==4
replace crc_fh_10yr_gp=f_fh5_med if sex==0 & crc_fh_10yr_gp==5
replace crc_fh_10yr_gp=f_fh6_med if sex==0 & crc_fh_10yr_gp==6

* current FH and PRS - women

tabstat crc_fh_prs_10yr if sex==0, by(crc_fh_prs_10yr_gp) ///
	statistics(p50) nototal save
matrix f_fhprs1=r(Stat1)
matrix f_fhprs2=r(Stat2)
matrix f_fhprs3=r(Stat3)
matrix f_fhprs4=r(Stat4)
matrix f_fhprs5=r(Stat5)
matrix f_fhprs6=r(Stat6)

scalar f_fhprs1_med=f_fhprs1[1,1]
scalar f_fhprs2_med=f_fhprs2[1,1]
scalar f_fhprs3_med=f_fhprs3[1,1]
scalar f_fhprs4_med=f_fhprs4[1,1]
scalar f_fhprs5_med=f_fhprs5[1,1]
scalar f_fhprs6_med=f_fhprs6[1,1]

replace crc_fh_prs_10yr_gp=f_fhprs1_med if sex==0 & crc_fh_prs_10yr_gp==1
replace crc_fh_prs_10yr_gp=f_fhprs2_med if sex==0 & crc_fh_prs_10yr_gp==2
replace crc_fh_prs_10yr_gp=f_fhprs3_med if sex==0 & crc_fh_prs_10yr_gp==3
replace crc_fh_prs_10yr_gp=f_fhprs4_med if sex==0 & crc_fh_prs_10yr_gp==4
replace crc_fh_prs_10yr_gp=f_fhprs5_med if sex==0 & crc_fh_prs_10yr_gp==5
replace crc_fh_prs_10yr_gp=f_fhprs6_med if sex==0 & crc_fh_prs_10yr_gp==6

* new FH and PRS - women

tabstat crc_newfh_prs_10yr if sex==0, by(crc_newfh_prs_10yr_gp) ///
	statistics(p50) nototal save
matrix f_newfhprs1=r(Stat1)
matrix f_newfhprs2=r(Stat2)
matrix f_newfhprs3=r(Stat3)
matrix f_newfhprs4=r(Stat4)
matrix f_newfhprs5=r(Stat5)
matrix f_newfhprs6=r(Stat6)

scalar f_newfhprs1_med=f_newfhprs1[1,1]
scalar f_newfhprs2_med=f_newfhprs2[1,1]
scalar f_newfhprs3_med=f_newfhprs3[1,1]
scalar f_newfhprs4_med=f_newfhprs4[1,1]
scalar f_newfhprs5_med=f_newfhprs5[1,1]
scalar f_newfhprs6_med=f_newfhprs6[1,1]

replace crc_newfh_prs_10yr_gp=f_newfhprs1_med if sex==0 & crc_newfh_prs_10yr_gp==1
replace crc_newfh_prs_10yr_gp=f_newfhprs2_med if sex==0 & crc_newfh_prs_10yr_gp==2
replace crc_newfh_prs_10yr_gp=f_newfhprs3_med if sex==0 & crc_newfh_prs_10yr_gp==3
replace crc_newfh_prs_10yr_gp=f_newfhprs4_med if sex==0 & crc_newfh_prs_10yr_gp==4
replace crc_newfh_prs_10yr_gp=f_newfhprs5_med if sex==0 & crc_newfh_prs_10yr_gp==5
replace crc_newfh_prs_10yr_gp=f_newfhprs6_med if sex==0 & crc_newfh_prs_10yr_gp==6

* new multivariable - women

tabstat crc_new_10yr if sex==0, by(crc_new_10yr_gp) ///
	statistics(p50) nototal save
matrix f_new1=r(Stat1)
matrix f_new2=r(Stat2)
matrix f_new3=r(Stat3)
matrix f_new4=r(Stat4)
matrix f_new5=r(Stat5)
matrix f_new6=r(Stat6)

scalar f_new1_med=f_new1[1,1]
scalar f_new2_med=f_new2[1,1]
scalar f_new3_med=f_new3[1,1]
scalar f_new4_med=f_new4[1,1]
scalar f_new5_med=f_new5[1,1]
scalar f_new6_med=f_new6[1,1]

replace crc_new_10yr_gp=f_new1_med if sex==0 & crc_new_10yr_gp==1
replace crc_new_10yr_gp=f_new2_med if sex==0 & crc_new_10yr_gp==2
replace crc_new_10yr_gp=f_new3_med if sex==0 & crc_new_10yr_gp==3
replace crc_new_10yr_gp=f_new4_med if sex==0 & crc_new_10yr_gp==4
replace crc_new_10yr_gp=f_new5_med if sex==0 & crc_new_10yr_gp==5
replace crc_new_10yr_gp=f_new6_med if sex==0 & crc_new_10yr_gp==6


*** MEN ***

* average - men

tabstat crc_av_10yr if sex==1, by(crc_av_10yr_gp) ///
	statistics(p50) nototal save
matrix m_av1=r(Stat1)
matrix m_av2=r(Stat2)
matrix m_av3=r(Stat3)
matrix m_av4=r(Stat4)
matrix m_av5=r(Stat5)
matrix m_av6=r(Stat6)

scalar m_av1_med=m_av1[1,1]
scalar m_av2_med=m_av2[1,1]
scalar m_av3_med=m_av3[1,1]
scalar m_av4_med=m_av4[1,1]
scalar m_av5_med=m_av5[1,1]
scalar m_av6_med=m_av6[1,1]

replace crc_av_10yr_gp=m_av1_med if sex==1 & crc_av_10yr_gp==1
replace crc_av_10yr_gp=m_av2_med if sex==1 & crc_av_10yr_gp==2
replace crc_av_10yr_gp=m_av3_med if sex==1 & crc_av_10yr_gp==3
replace crc_av_10yr_gp=m_av4_med if sex==1 & crc_av_10yr_gp==4
replace crc_av_10yr_gp=m_av5_med if sex==1 & crc_av_10yr_gp==5
replace crc_av_10yr_gp=m_av6_med if sex==1 & crc_av_10yr_gp==6

* current FH - men

tabstat crc_fh_10yr if sex==1, by(crc_fh_10yr_gp) ///
	statistics(p50) nototal save
matrix m_fh1=r(Stat1)
matrix m_fh2=r(Stat2)
matrix m_fh3=r(Stat3)
matrix m_fh4=r(Stat4)
matrix m_fh5=r(Stat5)
matrix m_fh6=r(Stat6)

scalar m_fh1_med=m_fh1[1,1]
scalar m_fh2_med=m_fh2[1,1]
scalar m_fh3_med=m_fh3[1,1]
scalar m_fh4_med=m_fh4[1,1]
scalar m_fh5_med=m_fh5[1,1]
scalar m_fh6_med=m_fh6[1,1]

replace crc_fh_10yr_gp=m_fh1_med if sex==1 & crc_fh_10yr_gp==1
replace crc_fh_10yr_gp=m_fh2_med if sex==1 & crc_fh_10yr_gp==2
replace crc_fh_10yr_gp=m_fh3_med if sex==1 & crc_fh_10yr_gp==3
replace crc_fh_10yr_gp=m_fh4_med if sex==1 & crc_fh_10yr_gp==4
replace crc_fh_10yr_gp=m_fh5_med if sex==1 & crc_fh_10yr_gp==5
replace crc_fh_10yr_gp=m_fh6_med if sex==1 & crc_fh_10yr_gp==6

* current FH and PRS - men

tabstat crc_fh_prs_10yr if sex==1, by(crc_fh_prs_10yr_gp) ///
	statistics(p50) nototal save
matrix m_fhprs1=r(Stat1)
matrix m_fhprs2=r(Stat2)
matrix m_fhprs3=r(Stat3)
matrix m_fhprs4=r(Stat4)
matrix m_fhprs5=r(Stat5)
matrix m_fhprs6=r(Stat6)

scalar m_fhprs1_med=m_fhprs1[1,1]
scalar m_fhprs2_med=m_fhprs2[1,1]
scalar m_fhprs3_med=m_fhprs3[1,1]
scalar m_fhprs4_med=m_fhprs4[1,1]
scalar m_fhprs5_med=m_fhprs5[1,1]
scalar m_fhprs6_med=m_fhprs6[1,1]

replace crc_fh_prs_10yr_gp=m_fhprs1_med if sex==1 & crc_fh_prs_10yr_gp==1
replace crc_fh_prs_10yr_gp=m_fhprs2_med if sex==1 & crc_fh_prs_10yr_gp==2
replace crc_fh_prs_10yr_gp=m_fhprs3_med if sex==1 & crc_fh_prs_10yr_gp==3
replace crc_fh_prs_10yr_gp=m_fhprs4_med if sex==1 & crc_fh_prs_10yr_gp==4
replace crc_fh_prs_10yr_gp=m_fhprs5_med if sex==1 & crc_fh_prs_10yr_gp==5
replace crc_fh_prs_10yr_gp=m_fhprs6_med if sex==1 & crc_fh_prs_10yr_gp==6

* new FH and PRS - men

tabstat crc_newfh_prs_10yr if sex==1, by(crc_newfh_prs_10yr_gp) ///
	statistics(p50) nototal save
matrix m_newfhprs1=r(Stat1)
matrix m_newfhprs2=r(Stat2)
matrix m_newfhprs3=r(Stat3)
matrix m_newfhprs4=r(Stat4)
matrix m_newfhprs5=r(Stat5)
matrix m_newfhprs6=r(Stat6)

scalar m_newfhprs1_med=m_newfhprs1[1,1]
scalar m_newfhprs2_med=m_newfhprs2[1,1]
scalar m_newfhprs3_med=m_newfhprs3[1,1]
scalar m_newfhprs4_med=m_newfhprs4[1,1]
scalar m_newfhprs5_med=m_newfhprs5[1,1]
scalar m_newfhprs6_med=m_newfhprs6[1,1]

replace crc_newfh_prs_10yr_gp=m_newfhprs1_med if sex==1 & crc_newfh_prs_10yr_gp==1
replace crc_newfh_prs_10yr_gp=m_newfhprs2_med if sex==1 & crc_newfh_prs_10yr_gp==2
replace crc_newfh_prs_10yr_gp=m_newfhprs3_med if sex==1 & crc_newfh_prs_10yr_gp==3
replace crc_newfh_prs_10yr_gp=m_newfhprs4_med if sex==1 & crc_newfh_prs_10yr_gp==4
replace crc_newfh_prs_10yr_gp=m_newfhprs5_med if sex==1 & crc_newfh_prs_10yr_gp==5
replace crc_newfh_prs_10yr_gp=m_newfhprs6_med if sex==1 & crc_newfh_prs_10yr_gp==6

* new multivariable - men

tabstat crc_new_10yr if sex==1, by(crc_new_10yr_gp) ///
	statistics(p50) nototal save
matrix m_new1=r(Stat1)
matrix m_new2=r(Stat2)
matrix m_new3=r(Stat3)
matrix m_new4=r(Stat4)
matrix m_new5=r(Stat5)
matrix m_new6=r(Stat6)

scalar m_new1_med=m_new1[1,1]
scalar m_new2_med=m_new2[1,1]
scalar m_new3_med=m_new3[1,1]
scalar m_new4_med=m_new4[1,1]
scalar m_new5_med=m_new5[1,1]
scalar m_new6_med=m_new6[1,1]

replace crc_new_10yr_gp=m_new1_med if sex==1 & crc_new_10yr_gp==1
replace crc_new_10yr_gp=m_new2_med if sex==1 & crc_new_10yr_gp==2
replace crc_new_10yr_gp=m_new3_med if sex==1 & crc_new_10yr_gp==3
replace crc_new_10yr_gp=m_new4_med if sex==1 & crc_new_10yr_gp==4
replace crc_new_10yr_gp=m_new5_med if sex==1 & crc_new_10yr_gp==5
replace crc_new_10yr_gp=m_new6_med if sex==1 & crc_new_10yr_gp==6

* some tables

tab age_gp5 crc_av_10yr_gp if sex==0
tab age_gp5 crc_fh_10yr_gp if sex==0
tab age_gp5 crc_fh_prs_10yr_gp if sex==0
tab age_gp5 crc_newfh_prs_10yr_gp if sex==0
tab age_gp5 crc_new_10yr_gp if sex==0

tab age_gp5 crc_av_10yr_gp if sex==1
tab age_gp5 crc_fh_10yr_gp if sex==1
tab age_gp5 crc_fh_prs_10yr_gp if sex==1
tab age_gp5 crc_newfh_prs_10yr_gp if sex==1
tab age_gp5 crc_new_10yr_gp if sex==1

* restricted to 10 years of follow-up
stset cr_sir_censor10_age, id(n_eid) enter(time age_calc) failure(cr10==1)

stsplit agegroup, at(40(5)80)

generate year = int(n_34_0_0 + int(_t0) + ((int(_t) - int(_t0))/2))
recode year (2005=2006) (2018/2019=2017)

merge m:1 agegroup sex year using "CR_incid_agegroup.dta"
drop if _merge==2 // ages outside the UKB age range
drop _merge


*** SIR analyses ***


*** compared to population rates ***

strate if sex==0, smr(cr_incid) 

strate age_gp10 if sex==0, smr(cr_incid) 

strate sex crc_av_10yr_gp if sex==0, smr(cr_incid) 
strate sex crc_fh_10yr_gp if sex==0, smr(cr_incid) 
strate sex crc_fh_prs_10yr_gp if sex==0, smr(cr_incid) 
strate sex crc_newfh_prs_10yr_gp if sex==0, smr(cr_incid) 
strate sex crc_new_10yr_gp if sex==0, smr(cr_incid) 


strate if sex==1, smr(cr_incid) 

strate age_gp10 if sex==1, smr(cr_incid) 


strate sex crc_av_10yr_gp if sex==1, smr(cr_incid) 
strate sex crc_fh_10yr_gp if sex==1, smr(cr_incid) 
strate sex crc_fh_prs_10yr_gp if sex==1, smr(cr_incid) 
strate sex crc_newfh_prs_10yr_gp if sex==1, smr(cr_incid) 
strate sex crc_new_10yr_gp if sex==1, smr(cr_incid) 

*** family history ***

strate sex cr_deg1 if sex==0, smr(cr_incid) 
strate sex cr_deg1 if sex==1, smr(cr_incid) 


*** graphs - women ***

colorpalette tab Hue circle
strate crc_av_10yr_gp if sex==0, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("Average risks: women", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:A}", size(*1.5) just(left)) ///
	name(sir_av_f, replace)

graph save "sir_av_f" "Graphs\SIR_av_f.gph", replace
graph export "Graphs\SIR_av_f.tif", name(sir_av_f) replace


colorpalette tab Hue circle
strate crc_fh_10yr_gp if sex==0, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("Family history alone model: women", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:B}", size(*1.5) just(left)) ///
	name(sir_fh_f, replace)

graph save "sir_fh_f" "Graphs\SIR_fh_f.gph", replace
graph export "Graphs\SIR_fh_f.tif", name(sir_fh_f) replace


colorpalette tab Hue circle
strate crc_fh_prs_10yr_gp if sex==0, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("Current family history and PRS model: women", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:C}", size(*1.5) just(left)) ///
	name(sir_fh_prs_f, replace)

graph save "sir_fh_prs_f" "Graphs\SIR_fh_prs_f.gph", replace
graph export "Graphs\SIR_fh_prs_f.tif", name(sir_fh_prs_f) replace



colorpalette tab Hue circle
strate crc_newfh_prs_10yr_gp if sex==0, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("New family history and PRS model: women", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:D}", size(*1.5) just(left)) ///
	name(sir_newfh_prs_f, replace)

graph save "sir_newfh_prs_f" "Graphs\SIR_newfh_prs_f.gph", replace
graph export "Graphs\SIR_newfh_prs_f.tif", name(sir_newfh_prs_f) replace

	
colorpalette tab Hue circle
strate crc_new_10yr_gp if sex==0, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("New multivariable model: women", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:E}", size(*1.5) just(left)) ///
	name(sir_new_f, replace)

graph save "sir_new_f" "Graphs\SIR_new_f.gph", replace
graph export "Graphs\SIR_new_f.tif", name(sir_new_f) replace	


* first-degree family history

colorpalette tab Hue circle
strate cr_deg1 if sex==0, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(-0.5 1.5)) xlabel(0 "No" 1 "Yes") ///
	xtitle("Affected first-degree relative", height(6)) ///
	title("First-degree family history: women", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.9 "{bf:F}", size(*1.5) just(left)) ///
	fxsize(40)  name(sir_deg1_f, replace)

graph save "sir_deg1_f" "Graphs\SIR_deg1_f.gph", replace
graph export "Graphs\SIR_deg1_f.tif", name(sir_deg1_f) replace	


* combined graph
graph combine "Graphs\SIR_av_f.gph" "Graphs\SIR_fh_f.gph"  ///
	"Graphs\SIR_fh_prs_f.gph" "Graphs\SIR_newfh_prs_f.gph" ///
	"Graphs\SIR_new_f.gph" "Graphs\SIR_deg1_f.gph", ///
	cols(2) iscale(0.5) xsize(7) ysize(9) graphregion(margin(small)) ///
	name(SIR_combined_f, replace)

graph save SIR_combined_f "Graphs\SIR_combined_f.gph", replace
graph export "Graphs\SIR_combined_f.tif", name(SIR_combined_f) replace	




*** graphs - men ***

* graphs - deciles 

colorpalette tab Hue circle
strate crc_av_10yr_gp if sex==1, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("Average risks: men", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:A}", size(*1.5) just(left)) ///
	name(sir_av_m, replace)

graph save "sir_av_m" "Graphs\SIR_av_m.gph", replace
graph export "Graphs\SIR_av_m.tif", name(sir_av_m) replace


colorpalette tab Hue circle
strate crc_fh_10yr_gp if sex==1, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("Family history alone model: men", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:B}", size(*1.5) just(left)) ///
	name(sir_fh_m, replace)

graph save "sir_fh_m" "Graphs\SIR_fh_m.gph", replace
graph export "Graphs\SIR_fh_m.tif", name(sir_fh_m) replace


colorpalette tab Hue circle
strate crc_fh_prs_10yr_gp if sex==1, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("Current family history and PRS model: men", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:C}", size(*1.5) just(left)) ///
	name(sir_fh_prs_m, replace)

graph save "sir_fh_prs_m" "Graphs\SIR_fh_prs_m.gph", replace
graph export "Graphs\SIR_fh_prs_m.tif", name(sir_fh_prs_m) replace

colorpalette tab Hue circle
strate crc_newfh_prs_10yr_gp if sex==1, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("New family history and PRS model: men", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:D}", size(*1.5) just(left)) ///
	name(sir_newfh_prs_m, replace)

graph save "sir_newfh_prs_m" "Graphs\SIR_newfh_prs_m.gph", replace
graph export "Graphs\SIR_newfh_prs_m.tif", name(sir_newfh_prs_m) replace

	
colorpalette tab Hue circle
strate crc_new_10yr_gp if sex==1, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(0 0.045)) xlabel(0 "0" 0.005 "0.5" 0.01 "1.0" 0.015 "1.5" ///
		0.02 "2.0" 0.025 "2.5" 0.03 "3.0" 0.035 "3.5" 0.04 "4.0" 0.045 "4.5") ///
	xtitle("10-year risk (%)", height(6)) ///
	title("New multivariable model: men", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.0065 "{bf:E}", size(*1.5) just(left)) ///
	name(sir_new_m, replace)

graph save "sir_new_m" "Graphs\SIR_new_m.gph", replace
graph export "Graphs\SIR_new_m.tif", name(sir_new_m) replace


* first-degree family history

colorpalette tab Hue circle
strate cr_deg1 if sex==1, smr(cr_incid) graph msymbol(O) ///
		mcolor("`r(p19)'") msize(large) ciopts(lcolor("`r(p19)'") lwidth(thick)) ///
	yline(1, lwidth(vthin) lcolor(gs8)) ///
	yscale(range(0.3 1.8)) ylabel(0.4(0.2)1.8, format(%02.1f) angle(h)) ///
	ytitle("Standardised incidence ratio", height(4)) ///
	xscale(range(-0.5 1.5)) xlabel(0 "No" 1 "Yes") ///
	xtitle("Affected first-degree relative", height(6)) ///
	title("First-degree family history: men", size(*0.75)) ///
	legend(on cols(1) ring(0) bplacement(nwest) region(lwidth(none)) ///
		order(2 1) label(1 "95% CI") label(2 "SIR") size(*0.80) symxsize(7) ///
		keygap(2) rowgap(*0.3)) ///
	text(1.95 -0.9 "{bf:F}", size(*1.5) just(left)) ///
	fxsize(40)  name(sir_deg1_m, replace)

graph save "sir_deg1_m" "Graphs\SIR_deg1_m.gph", replace
graph export "Graphs\SIR_deg1_m.tif", name(sir_deg1_m) replace	


* combined graph
graph combine "Graphs\SIR_av_m.gph" "Graphs\SIR_fh_m.gph"  ///
	"Graphs\SIR_fh_prs_m.gph" "Graphs\SIR_newfh_prs_m.gph" ///
	"Graphs\SIR_new_m.gph" "Graphs\SIR_deg1_m.gph", ///
	cols(2) iscale(0.5) xsize(7) ysize(9) ///
	graphregion(margin(small)) name(SIR_combined_m, replace)

graph save SIR_combined_m "Graphs\SIR_combined_m.gph", replace
graph export "Graphs\SIR_combined_m.tif", name(SIR_combined_m) replace		


* combined graph JAMA
graph combine "Graphs\SIR_newfh_prs_f.gph" "Graphs\SIR_new_f.gph" ///
	"Graphs\SIR_deg1_f.gph" "Graphs\SIR_newfh_prs_m.gph" ///
	"Graphs\SIR_new_m.gph" "Graphs\SIR_deg1_m.gph", cols(3) iscale(0.5) ///
	xsize(7) graphregion(margin(small)) name(SIR_combined_jama, replace)
graph save "SIR_combined_jama" "Graphs\SIR_combined_jama.gph", replace
graph export "Graphs\SIR_combined_jama.tif", name(SIR_combined_jama) replace	


log close sir
