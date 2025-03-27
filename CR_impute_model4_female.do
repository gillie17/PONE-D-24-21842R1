*** UK Biobank -- colorectal cancer modelling using imputed data ***

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Model4_female_imputed_`filedate'.log", ///
	name(model4f_impute) replace

use "CR_analysis_imputed.dta", clear

mi convert flong // do this if doing separate analyses on the _mi_m

*** MODEL 4 - Female ***

**# base model ***

mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi_c ///
	ln_activity_mets i.alcohol if sex==0
mi estimate, coeflegend 
	
* to do the contrast = F test is the overall association for that variable
mi test 1.alcohol 2.alcohol 3.alcohol

* drop bmi
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ///
	ln_activity_mets i.alcohol if sex==0
mi test 1.alcohol 2.alcohol 3.alcohol

* drop diabetes
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever ln_activity_mets ///
	i.alcohol if sex==0
mi test 1.alcohol 2.alcohol 3.alcohol

* drop alcohol
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever ln_activity_mets if sex==0		

* drop activity
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever if sex==0	
	

* final base model
mi estimate, hr saving(models\model4f_base, replace): stcox std_prs_xb ///
	cr_deg1 smoke_ever if sex==0	
	

**# add diet variables to the model ***

* forwards selection - add variables in same order as for men

* dried fruit; processed meat; cereal; beef; pork; raw veg; fresh fruit;
* white bread; whole bread; cooked veg

* dried fruit - yes/no - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever ///
	diet_fruit_dried_cat if sex==0	

* processed meat - use categories because of answer format - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever i.diet_procmeat_cat ///
	if sex==0	
mi test 1.diet_procmeat_cat 2.diet_procmeat_cat 3.diet_procmeat_cat

* cereal - use categories because of distribution - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever i.diet_cereal_cat if sex==0	
mi test 1.diet_cereal_cat 2.diet_cereal_cat 3.diet_cereal_cat

* beef - use categories because of answer format - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever i.diet_beef_cat if sex==0	
mi test 1.diet_beef_cat 2.diet_beef_cat 3.diet_beef_cat
	
* pork - use categories because of answer format - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever i.diet_pork_cat if sex==0	
mi test 1.diet_pork_cat 2.diet_pork_cat 3.diet_pork_cat	

* raw vegetables - use continuous measure - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever diet_veg_raw if sex==0	

* fresh fruit - use continuous measure - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever diet_fruit_fresh if sex==0	
	
* white bread - use categories because of distribution - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever ///
	i.diet_bread_white_cat_new if sex==0	
mi test 1.diet_bread_white_cat_new 2.diet_bread_white_cat_new ///
	3.diet_bread_white_cat_new

* whole bread - use categories because of distribution - drop 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever ///
	i.diet_bread_whole_cat_new if sex==0	
mi test 1.diet_bread_whole_cat_new 2.diet_bread_whole_cat_new ///
	3.diet_bread_whole_cat_new

* cooked vegetables - use continuous measure - drop	
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever diet_veg_cook if sex==0	


* model with diet risk factors

mi estimate, hr saving(models\model4f_diet, replace): stcox std_prs_xb ///
	cr_deg1 smoke_ever if sex==0


*** add other risk factors to the model ***
	
* blood markers together 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30760_0_0_c ///
	n_30870_0_0_c n_30780_0_0_c n_30690_0_0_c if sex==0
	
* drop n_30780_0_0 
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30760_0_0_c ///
	n_30870_0_0_c n_30690_0_0_c if sex==0
	
* drop n_30760_0_0
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	n_30690_0_0_c if sex==0

* drop n_30690_0_0
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c if sex==0	

* NSAID - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c nsaid if sex==0
	
* calcium - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	supp_calcium if sex==0

* vitamin D - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	supp_vitamin_d if sex==0

* fishoil - yes/no - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c fishoil ///
	if sex==0
	
* screening - keep yes/no but drop time
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever ///
	n_30870_0_0_c screen10_new screen10_new_time if sex==0
	
* menopause/HRT - drop
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	screen10_new i.menop_hrt if sex==0
mi test 1.menop_hrt 2.menop_hrt


*** selected model ***	
	
mi estimate, hr saving(models\model4f_selected, replace): stcox std_prs_xb ///
	cr_deg1 smoke_ever n_30870_0_0_c screen10_new if sex==0


*** check for time-varying risk factors = proportional hazards test ***

* this is the same as doing a proportional hazards test
mi estimate, hr: stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c ///
	screen10_new if sex==0, tvc(std_prs_xb cr_deg1 smoke_ever n_30870_0_0 ///
	screen10_new)	

* proportional hazards test in first imputation dataset - no problems

stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c screen10_new ///
	if sex==0 & _mi_m==1, hr
estat phtest, detail

*** FINAL MODEL 4 - female ***

mi estimate, saving(models\model4f_final, replace): stcox std_prs_xb ///
	cr_deg1 smoke_ever n_30870_0_0_c screen10_new if sex==0
mi estimate, hr base
mi estimate, coeflegend 

how_many_imputations, cv_se(0.05) 

*** model fit ***

mi extract 1

keep if sex==0

stcox std_prs_xb cr_deg1 smoke_ever n_30870_0_0_c screen10_new, nohr mgale(mg)

predict cs, csnell // Cox-Snell residuals

drop mg

stset cs, failure(cr==1)
sts generate H = na // Nelson-Aalen cumulative hazard
twoway line H cs if cs<=0.05, sort lcolor(orange) || ///
	line cs cs if cs<=0.05, sort lcolor(ebblue) ///
	aspect(1) xsize(3.5) xtitle(, margin(t +2)) ///
	xscale(range(0 0.05)) xlabels(0(0.01)0.05, format(%03.2f)) ///
	yscale(range(0 0.05)) ylabels(0(0.01)0.05, format(%03.2f) angle(h))  ///
	legend(cols(1) size(small) region(lwidth(none))) ///
		name(NelsonAalen_f, replace)
graph save "NelsonAalen_f" "Graphs\CoxSnell_NelsonAalen_f.gph", replace
graph export "Graphs\CoxSnell_NelsonAalen_f.tif", as(tif) ///
	name("NelsonAalen_f") replace


log close model4f_impute


