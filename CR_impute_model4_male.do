*** UK Biobank -- colorectal cancer modelling using imputed data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Model4_male_imputed_`filedate'.log", ///
	name(model4m_impute) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m


*** MODEL 4 - male ***


**# base model ***

mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	ln_activity_mets i.alcohol if sex==1
mi estimate, coeflegend 
	
* to do the contrast = F test is the overall association for that variable
mi test 1.alcohol 2.alcohol 3.alcohol

* drop ln_activity_mets - this is the final base model
mi estimate, hr saving(models\model4m_base, replace): stcox std_prs_xb ///
	cr_deg1 diabetes smoke_ever ln_bmi_c i.alcohol if sex==1

* check the linear combinations
mi estimate (alc12: _b[1.alcohol] - _b[2.alcohol]) ///
	(alc13: _b[1.alcohol] - _b[3.alcohol]) ///
	(alc23: _b[2.alcohol] - _b[3.alcohol]) using models\model4m_base	
	

**# add diet variables to the model ***

* forwards selection - add variables in order of univariate association 

* dried fruit; processed meat; cereal; beef; pork; fishoil; raw veg; fresh fruit
* white bread; whole bread; cooked veg 
 
* dried fruit - yes/no - keep
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol diet_fruit_dried_cat if sex==1

* processed meat - use categories because of answer format - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol diet_fruit_dried_cat i.diet_procmeat_cat if sex==1
mi test 1.diet_procmeat_cat 2.diet_procmeat_cat 3.diet_procmeat_cat
	
* cereal - use categories because of distribution - keep but drop dried fruit
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol diet_fruit_dried_cat i.diet_cereal_cat if sex==1
mi test 1.diet_cereal_cat 2.diet_cereal_cat 3.diet_cereal_cat
	
* beef - use categories because of answer format - keep 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat if sex==1
mi test 1.diet_beef_cat 2.diet_beef_cat 3.diet_beef_cat
	
* pork - use categories because of answer format - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_pork_cat if sex==1
mi test 1.diet_pork_cat 2.diet_pork_cat 3.diet_pork_cat	
mi test 1.diet_beef_cat 2.diet_beef_cat 3.diet_beef_cat

* raw vegetables - use continuous measure - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat diet_veg_raw if sex==1

* fresh fruit - use continuous measure - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat diet_fruit_fresh if sex==1

* white bread - use categories because of distribution - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_bread_white_cat_new ///
	if sex==1
mi test 1.diet_bread_white_cat_new 2.diet_bread_white_cat_new ///
	3.diet_bread_white_cat_new

* whole bread - use categories because of distribution - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_bread_whole_cat_new ///
	if sex==1
mi test 1.diet_bread_whole_cat_new 2.diet_bread_whole_cat_new ///
	3.diet_bread_whole_cat_new

* cooked vegetables - use continuous measure - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat diet_veg_cook if sex==1

**# check the dropped variables

* dried fruit - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol  i.diet_cereal_cat i.diet_beef_cat diet_fruit_dried_cat if sex==1

* processed meat - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_procmeat_cat if sex==1

* pork - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_pork_cat if sex==1

* raw veg - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat diet_veg_raw if sex==1

* fresh fruit - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat diet_fruit_fresh if sex==1

* white bread - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_bread_white_cat_new ///
	if sex==1	
	
* whole bread - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.diet_bread_whole_cat_new ///
	if sex==1

* cooked veg - don't add
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat diet_veg_cook if sex==1
	
* model with diet risk factors

mi estimate, hr saving(models\model4m_diet, replace): stcox std_prs_xb ///
	cr_deg1 diabetes smoke_ever ln_bmi_c i.alcohol i.diet_cereal_cat ///
	i.diet_beef_cat if sex==1

**# add other risk factors to the model ***

* blood markers - drop n_30780_0_0 
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat n_30760_0_0 n_30870_0_0 ///
	n_30780_0_0 n_30690_0_0 if sex==1

* drop n_30760_0_0
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat n_30760_0_0 n_30870_0_0 ///
	n_30690_0_0 if sex==1

* drop n_30690_0_0
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat n_30870_0_0 n_30690_0_0 ///
	if sex==1

* drop n_30690_0_0
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat n_30870_0_0 if sex==1

* NSAID - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat nsaid if sex==1
	
* calcium - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat supp_calcium if sex==1

* vitamin D - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat supp_vitamin_d if sex==1	

* fishoil - yes/no - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat fishoil if sex==1
	
* screening	- keep yes/no and the time variable
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat screen10_new screen10_new_time ///
	if sex==1	
	
*** selected model *** 

mi estimate, hr saving(models\model4m_selected, replace): stcox std_prs_xb ///
	cr_deg1 diabetes smoke_ever ln_bmi_c i.alcohol i.diet_cereal_cat ///
	i.diet_beef_cat screen10_new screen10_new_time if sex==1

*** check for time-varying risk factors  = proportional hazards test ***

* this is the same as doing a proportional hazards test
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat screen10_new ///
	screen10_new_time if sex==1, tvc(std_prs_xb cr_deg1 diabetes smoke_ever ///
	ln_bmi_c i.alcohol i.diet_cereal_cat i.diet_beef_cat screen10_new ///
	screen10_new_time)

* proportional hazards test in first imputation dataset
* none of these graphs show a clear problem; the test is very sensitive	
stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c i.alcohol ///
	i.diet_cereal_cat i.diet_beef_cat screen10_new screen10_new_time ///
	if sex==1 & _mi_m==1
estat phtest, detail
estat phtest, plot(std_prs_xb) name(std_prs_xb, replace)
estat phtest, plot(cr_deg1) name(cr_deg1, replace)
estat phtest, plot(diabetes) name(diabetes, replace)
estat phtest, plot(1.alcohol) name(alcohol1, replace)
estat phtest, plot(2.alcohol) name(alcohol2, replace)
estat phtest, plot(3.alcohol) name(alcohol3, replace)
estat phtest, plot(screen10_new) name(screen10_new, replace)
estat phtest, plot(screen10_new_time) name(screen10_new_time, replace)
	
*** FINAL MODEL 4 - male ***

mi estimate, saving(models\model4m_final, replace): stcox std_prs_xb ///
	cr_deg1 diabetes smoke_ever ln_bmi_c i.alcohol i.diet_cereal_cat ///
	i.diet_beef_cat screen10_new screen10_new_time if sex==1
mi estimate, hr base
mi estimate, coeflegend 

how_many_imputations, cv_se(0.05) 

*** model fit ***

mi extract 1

keep if sex==1

stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c i.alcohol ///
	i.diet_cereal_cat i.diet_beef_cat screen10_new screen10_new_time, ///
	nohr mgale(mg)

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
		name(NelsonAalen_m, replace)
graph save "NelsonAalen_m" "Graphs\CoxSnell_NelsonAalen_m.gph", replace
graph export "Graphs\CoxSnell_NelsonAalen_m.tif", as(tif) ///
	name("NelsonAalen_m") replace


log close model4m_impute



