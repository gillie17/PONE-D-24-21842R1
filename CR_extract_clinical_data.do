*** UK Biobank -- colorectal cancer dataset extraction -- used by CR_read_data.do ***

version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\UKB Data"

local usevars n_eid /// ID number
n_20002_0_0-n_20002_0_28 /// non-cancer illness code 
n_20008_0_0-n_20008_0_28 /// non-cancer illness year
n_20009_0_0-n_20009_0_28 /// non-cancer illness age
n_20004_0_0-n_20004_0_31 /// operation code
n_20010_0_0-n_20010_0_31 /// operation code year
n_20011_0_0-n_20011_0_31 /// operation code age
s_41271_0_0-s_41271_0_46 /// ICD9 code
ts_41281_0_0-ts_41281_0_46 /// ICD9 date
s_41270_0_0-s_41270_0_212 /// ICD10 code
ts_41280_0_0-ts_41280_0_212 /// ICD10 date
s_41272_0_0-s_41272_0_116 /// OPSC4 operation code
ts_41282_0_0-ts_41282_0_116 /// OPSC4 operation date
n_2345_0_0 n_2355_0_0  /// screening  questions
n_50_0_0 /// height
n_21002_0_0 /// weight
n_21001_0_0 /// BMI
n_2724_0_0 /// menopause
n_2814_0_0 /// HRT ever
n_3536_0_0 /// age started HRT
n_3546_0_0 /// age last used HRT
n_20116_0_0 /// smoking status
n_20161_0_0 /// pack years of smoking
n_2897_0_0 /// age stopped smoking
n_1558_0_0 /// drinking status
n_2443_0_0 /// diabetes diagnosed by dr
n_2976_0_0 // age diabetes diagnosed

use `usevars' using "UKB_ClinicalVariables.dta", clear

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"
save "CR_clinical_1.dta", replace


cd "R:\Reference\UKBPhenotype\Basket 43105 (2008311)\Stata"

local usevars n_eid /// ID number
n_6154_0_0-n_6154_0_5 /// medication
n_6155_0_0-n_6155_0_6 /// vitamin and mineral supplements
n_6179_0_0-n_6179_0_5 /// mineral and other dietarty supplements
n_864_0_0 /// walk times per week
n_874_0_0 /// duration of walks
n_884_0_0 /// moderate activity times per week
n_894_0_0 /// duration of moderate activity
n_904_0_0 /// vigorous activity times per week
n_914_0_0 /// duration of vigorous activity
n_924_0_0 /// usual walking pace
n_20003_0_0-n_20003_0_47 /// medications
n_1289_0_0 /// Cooked vegetable intake
n_1299_0_0 /// Salad/raw vegetable intake
n_1309_0_0 /// Fresh fruit intake
n_1319_0_0 /// Dried fruit intake
n_1329_0_0 /// Oily fish intake
n_1339_0_0 /// Non-oily fish intake
n_1349_0_0 /// Processed meat intake
n_1359_0_0 /// Poultry intake
n_1369_0_0 /// Beef intake
n_1379_0_0 /// Lamb/mutton intake
n_1389_0_0 /// Pork intake
*n_3680_0_0 /// Age when last ate meat
*n_6144_0_0 /// Never eat eggs, dairy, wheat, sugar
n_1408_0_0 /// Cheese intake
*n_1418_0_0 /// Milk type used
*n_1428_0_0 /// Spread type
*n_2654_0_0 /// Non-butter spread type details
n_1438_0_0 /// Bread intake
*n_1448_0_0 /// Bread type
n_1458_0_0 /// Cereal intake
*n_1468_0_0 /// Cereal type
*n_1478_0_0 /// Salt added to food
*n_1488_0_0 /// Tea intake
*n_1508_0_0 /// Coffee type
*n_1518_0_0 /// Hot drink temperature
n_1528_0_0 /// Water intake
n_1538_0_0 /// Major dietary changes in the last 5 years
n_1548_0_0 // Variation in diet

use `usevars' using "ClinicalData_43105.dta", clear

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"
save "CR_clinical_2.dta", replace


cd "R:\Reference\UKBPhenotype\Basket 43105 (2008311)\Stata"

local usevars n_eid /// ID number
n_6154_0_0-n_6154_0_5 /// medication
n_6155_0_0-n_6155_0_6 /// vitamin and mineral supplements
n_6179_0_0-n_6179_0_5 /// mineral and other dietarty supplements
n_884_0_0 /// moderate excercise times per week
n_904_0_0 /// vigorous excercise times per week
n_924_0_0 /// usual walking pace
n_6164_0_0-n_6164_0_4 /// types of physical activity in last 4 weeks
n_20003_0_0-n_20003_0_47 /// medications
n_1289_0_0 /// Cooked vegetable intake
n_1299_0_0 /// Salad/raw vegetable intake
n_1309_0_0 /// Fresh fruit intake
n_1319_0_0 /// Dried fruit intake
n_1329_0_0 /// Oily fish intake
n_1339_0_0 /// Non-oily fish intake
n_1349_0_0 /// Processed meat intake
n_1359_0_0 /// Poultry intake
n_1369_0_0 /// Beef intake
n_1379_0_0 /// Lamb/mutton intake
n_1389_0_0 /// Pork intake
*n_3680_0_0 /// Age when last ate meat
*n_6144_0_0 /// Never eat eggs, dairy, wheat, sugar
n_1408_0_0 /// Cheese intake
*n_1418_0_0 /// Milk type used
*n_1428_0_0 /// Spread type
*n_2654_0_0 /// Non-butter spread type details
n_1438_0_0 /// Bread intake
*n_1448_0_0 /// Bread type
n_1458_0_0 /// Cereal intake
*n_1468_0_0 /// Cereal type
*n_1478_0_0 /// Salt added to food
*n_1488_0_0 /// Tea intake
*n_1508_0_0 /// Coffee type
*n_1518_0_0 /// Hot drink temperature
n_1528_0_0 /// Water intake
n_1538_0_0 /// Major dietary changes in the last 5 years
n_1548_0_0 // Variation in diet

use `usevars' using "ClinicalData_43105.dta", clear

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"
save "CR_clinical_2.dta", replace


cd "R:\Reference\UKBPhenotype\Basket 670395 (2014059)\Stata"

local usevars n_eid /// ID number
ts_131626_0_0 /// date K50 first reported (crohn's disease) 
ts_131628_0_0 /// date K51 first reported (ulcerative colitis) 
ts_131630_0_0 /// date K52 first reported (other ibd)
ts_130706_0_0 /// date E10 first reported (type 1 diabetes)
ts_130708_0_0 /// date E11 first reported (type 2 diabetes)
ts_130712_0_0 /// date E12 first reported (other specified diabetes)
ts_130714_0_0 // date E14 first reported (unspecified diabetes)

use `usevars' using "ClinicalData_ukb670395.dta", clear

format ts_131626_0_0 ts_131628_0_0 ts_131630_0_0 ts_130706_0_0 ts_130714_0_0 %td

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data"
save "CR_clinical_3.dta", replace

 





