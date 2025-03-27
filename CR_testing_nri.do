cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"

local filedate: display %tdCCYYNNDD td(`c(current_date)')
log using "Results\Testing_NRI_`filedate'.log", ///
	name(nri) replace

use "CR_testing_analysis.dta", clear

keep if ancestry_code==9

egen fh_1 = cut(crc_fh_10yr), at(0,0.01,1)
egen new_fhprs_1 = cut(crc_newfh_prs_10yr), at(0,0.01,1)
egen new_1 = cut(crc_new_10yr), at(0,0.01,1)

egen fh_2 = cut(crc_fh_10yr), at(0,0.02,1)
egen new_fhprs_2 = cut(crc_newfh_prs_10yr), at(0,0.02,1)
egen new_2 = cut(crc_new_10yr), at(0,0.02,1)

* 1% threshold
bysort cr10: tab cr_deg1 new_fhprs_1 if sex==0
bysort cr10: tab cr_deg1 new_fhprs_1 if sex==1

bysort cr10: tab cr_deg1 new_1 if sex==0
bysort cr10: tab cr_deg1 new_1 if sex==1

bysort cr10: tab fh_1 new_fhprs_1 if sex==0
bysort cr10: tab fh_1 new_fhprs_1 if sex==1

bysort cr10: tab fh_1 new_1 if sex==0
bysort cr10: tab fh_1 new_1 if sex==1

* 2% threshold
bysort cr10: tab cr_deg1 new_fhprs_2 if sex==0
bysort cr10: tab cr_deg1 new_fhprs_2 if sex==1

bysort cr10: tab cr_deg1 new_2 if sex==0
bysort cr10: tab cr_deg1 new_2 if sex==1

bysort cr10: tab fh_2 new_fhprs_2 if sex==0
bysort cr10: tab fh_2 new_fhprs_2 if sex==1

bysort cr10: tab fh_2 new_2 if sex==0
bysort cr10: tab fh_2 new_2 if sex==1


* three groups

egen fh_1_2 = cut(crc_fh_10yr), at(0,0.01,0.02,1)
egen new_fhprs_1_2 = cut(crc_newfh_prs_10yr), at(0,0.01,0.02,1)
egen new_1_2 = cut(crc_new_10yr), at(0,0.01,0.02,1)

bysort cr10: tab fh_1_2 new_fhprs_1_2 if sex==0
bysort cr10: tab fh_1_2 new_fhprs_1_2 if sex==1

bysort cr10: tab fh_1_2 new_1_2 if sex==0
bysort cr10: tab fh_1_2 new_1_2 if sex==1



log close nri