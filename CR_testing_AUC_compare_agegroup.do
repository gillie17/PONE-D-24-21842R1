
version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Testing_AUC_compare_agegroup_`filedate'.log", ///
	name(auc_compare) replace

use "CR_testing_analysis.dta", clear

keep if ancestry_code==9

* restricted to 10 years of follow-up
stset cr_cox_censor10_age, id(n_eid) enter(time age_calc) failure(cr10==1)


*** AUCs by 10-year age group

* women
bysort age_gp10: somersd cr10 crc_av_10yr crc_fh_10yr crc_fh_prs_10yr ///
	crc_newfh_prs_10yr crc_new_10yr if sex==0, transf(c)

* men
bysort age_gp10: somersd cr10 crc_av_10yr crc_fh_10yr crc_fh_prs_10yr ///
	crc_newfh_prs_10yr crc_new_10yr if sex==1, transf(c)

* compare within age groups	
	
foreach n in 1 2 3 {
foreach s in 0 1 {

di ""
di "age group = " "`n'"
di "sex = " "`s'"
	
somersd cr10 crc_av_10yr crc_fh_10yr crc_fh_prs_10yr crc_newfh_prs_10yr ///
	crc_new_10yr if sex==`s' & age_gp10==`n', transf(c)
lincom crc_newfh_prs_10yr-crc_new_10yr 
	
lincom crc_newfh_prs_10yr-crc_fh_prs_10yr 
lincom crc_newfh_prs_10yr-crc_fh_10yr 
lincom crc_newfh_prs_10yr-crc_av_10yr 

lincom crc_new_10yr-crc_fh_prs_10yr 
lincom crc_new_10yr-crc_fh_10yr 
lincom crc_new_10yr-crc_av_10yr 
}
}
	
*** program to compare AUCs between age groups

* arguments needed to run program are variable name and 0 or 1 (female and male)
* e.g. compare crc_new_10yr 0	
	
capture program drop compareauc
program compareauc

di _newline
di "`1'"
di "Sex = " "`2'" " (0=female, 1=male)"	
di _newline
di "40-49: test AUC=0.5"

quietly somersd cr10 `1' if sex==`2' & age_gp10==1, transf(c)
test `1'=0.5
scalar auc40=e(b)[1,1]
scalar var40=e(V)[1,1]

di _newline
di "50-59: test AUC=0.5"

quietly somersd cr10 `1' if sex==`2' & age_gp10==2, transf(c)

test `1'=0.5
scalar auc50=e(b)[1,1]
scalar var50=e(V)[1,1]

di _newline
di "60-69: test AUC=0.5"	
quietly somersd cr10 `1' if sex==`2' & age_gp10==3, transf(c)

test `1'=0.5
scalar auc60=e(b)[1,1]
scalar var60=e(V)[1,1]	
	
scalar z_40_50 = abs((auc50 - auc40) / sqrt(var40 + var50))
scalar z_50_60 = abs((auc50 - auc60) / sqrt(var50 + var60))
scalar z_40_60 = abs((auc40 - auc60) / sqrt(var40 + var60))

di _newline
di "Comparisons between age groups"
di _newline
di "40 vs 50:   z = " %05.3g z_40_50 "    p = "%05.3g 2*(1 - normal(z_40_50)) 
di "50 vs 60:   z = " %05.3g z_50_60 "    p = "%05.3g 2*(1 - normal(z_50_60)) 
di "40 vs 60:   z = " %05.3g z_40_60 "    p = "%05.3g 2*(1 - normal(z_40_60)) 

end


foreach v of varlist crc_av_10yr crc_fh_10yr crc_fh_prs_10yr ///
	crc_newfh_prs_10yr crc_new_10yr {
foreach n of numlist 0/1 {

compareauc `v' `n'

}
}

log close auc_compare

exit