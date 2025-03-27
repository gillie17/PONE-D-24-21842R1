
*** WOMEN n=149928)  ***

* menop/hrt dummy variables
power cox, alpha(0.05) n(91005) hratio(1.2) sd(0.48) eventprob(0.009)
power cox, alpha(0.05) n(95757) hratio(1.07) sd(0.48) eventprob(0.009)

* ln_bmi
power cox, alpha(0.05) n(152088) hratio(1.98) sd(0.18) eventprob(0.009)

* fruit - dried
power cox, alpha(0.05) n(150000) hratio(0.96) sd(0.5) eventprob(0.009)

* polyp
power cox, alpha(0.05) n(152088) hratio(1.34) sd(0.13) eventprob(0.009)


* blood markers
power cox, alpha(0.05) n(150000) hratio(0.93) sd(0.38) eventprob(0.009)
power cox, alpha(0.05) n(150000) hratio(1.06) sd(0.85) eventprob(0.009)
power cox, alpha(0.05) n(150000) hratio(0.90) sd(0.87) eventprob(0.009)
power cox, alpha(0.05) n(150000) hratio(1.11) sd(1.12) eventprob(0.009)

* diabetes
power cox, alpha(0.05) n(152008) hratio(1.25) sd(0.18) eventprob(0.009)


*** MEN n=127323 ***

* nsaid - KEEP

mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.screen15_cat ///
	nsaid if sex==1

tab nsaid cr if sex==1 &_mi_m==0, co

summ nsaid if sex==1 &_mi_m==0

power cox, alpha(0.05) n(127323) hratio(0.908) sd(0.47) eventprob(0.0143) 

power cox, alpha(0.05) n(127323) power(0.8) sd(0.47) eventprob(0.0143) ///
	effect(hr) direction(lower)

	
* dried fruit - ???
	
mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.screen15_cat ///
	diet_fruit_dried_cat if sex==1

summ diet_fruit_dried_cat if sex==1 &_mi_m==0	

power cox, alpha(0.05) n(127323) hratio(0.938) sd(0.46) eventprob(0.0143) 

power cox, alpha(0.05) n(127323) power(0.8) sd(0.46) eventprob(0.0143) ///
	effect(hr) direction(lower)




* pork

mi estimate, hr: stcox std_prs_xb cr_deg1 diabetes smoke_ever ln_bmi ///
	i.alcohol i.diet_cereal_cat i.diet_beef_cat i.screen15_cat ///
	i.diet_pork_cat if sex==1

tab diet_pork_cat, generate (diet_pork_cat_d)

summ diet_pork_cat_d1 diet_pork_cat_d2 diet_pork_cat_d3 if sex==1

* power
power cox, alpha(0.05) n(89968) hratio(0.936) sd(0.49) eventprob(0.0143) /// 
	effect(hr) 
power cox, alpha(0.05) n(46818) hratio(1.064) sd(0.43) eventprob(0.0143) /// 
	effect(hr) 
power cox, alpha(0.05) n(20259) hratio(1.179) sd(0.20) eventprob(0.0143) /// 
	effect(hr) 

* minimum effect size
power cox, alpha(0.05) n(89968) power(0.8) sd(0.49) eventprob(0.0143) ///
	effect(hr) direction(lower)
power cox, alpha(0.05) n(46818) power(0.8) sd(0.43) eventprob(0.0143) ///
	effect(hr) direction(upper)
power cox, alpha(0.05) n(20259) power(0.8) sd(0.20) eventprob(0.0143) ///
	effect(hr) direction(upper)




power cox, alpha(0.05) n(130064) hratio(1.083) sd(0.064) eventprob(0.0143)

* diabetes
power cox, alpha(0.05) n(130064) hratio(1.25) sd(0.24) eventprob(0.0143)



