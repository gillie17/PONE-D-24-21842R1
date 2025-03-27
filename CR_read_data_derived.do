*** UK Biobank -- dereived variable creation -- called by CR_read_data.do ***

generate fup5_time=fup_time
replace fup5_time=5 if fup_time>5

generate fup10_time=fup_time
replace fup10_time=10 if fup_time>10

generate age_calc_2=age_calc^2
label variable age_calc_2 "Age at baseline squared"

*** PRS

summ prs_rr
scalar sd_prs_rr=r(sd)
generate prs_rr_sd=prs_rr/sd_prs_rr
label variable prs_rr_sd "Per SD of PRS"

summ prs_xb
scalar sd_prs_xb=r(sd)
generate prs_xb_sd=prs_xb/sd_prs_xb
label variable prs_xb_sd "Per SD of log PRS"

* standardised PRS from lc_prs_xb

egen std_prs_xb=std(lc_prs_xb)
label variable std_prs_xb "Standardised linear combination PRS"


*** family history

replace cr_deg1n=. if miss_mum==1 & miss_dad==1 & miss_sib==1

generate cr_deg1=cr_deg1n
recode cr_deg1 (2 3=1)

label variable cr_deg1 "Any affected first-degree relative"
label values cr_deg1 yesno

generate cr_deg1_012=cr_deg1n
recode cr_deg1_012 (3=2)
label variable cr_deg1_012 "Affected first-degree relative"
label define fh_label 0 "None" 1 "One" 2 "Two or more"
label values cr_deg1_012 fh_label

*** body composition

generate height=n_50_0_0 
label variable height "Height (cm)"
generate weight=n_21002_0_0 
label variable weight "Weight (kg)"
generate bmi=n_21001_0_0
label variable bmi "BMI (kg/m^2)" 
generate ln_bmi=ln(bmi)
label variable ln_bmi "BMI (natural log)"


/* POLYPS BEFORE BASELINE */

*** self-reported rectal or colon adenoma/polyp

generate polyp_s=0
generate polyp_s_year=.
generate polyp_s_age=.

forvalues n=0/28 {
	
quietly replace polyp_s=1 if n_20002_0_`n'==1460
quietly replace polyp_s_year=n_20008_0_`n' if n_20002_0_`n'==1460 & ///
	n_20008_0_`n'<polyp_s_year & n_20008_0_`n'!=-1 & n_20008_0_`n'!=-3
quietly replace polyp_s_age=n_20009_0_`n' if n_20002_0_`n'==1460 & ///
	n_20009_0_`n'<polyp_s_age & n_20009_0_`n'!=-1 & n_20009_0_`n'!=-3
}

generate polyp_s_time=age_calc-polyp_s_age if polyp_s==1 
replace polyp_s_time=0 if polyp_s_time<0

label variable polyp_s "Self-reported polyp"
label values polyp_s yesno
label variable polyp_s_year "Self-reported polyp: year"
label variable polyp_s_age "Self-reported polyp: age"
label variable polyp_s_time "Self-reported polyp: time since"

*** self-reported polyp procedure

generate polyp_so=0
generate polyp_so_year=.
generate polyp_so_age=.

forvalues n=0/31 {
	
quietly replace polyp_so=1 if n_20004_0_`n'==1463
quietly replace polyp_so_year=n_20010_0_`n' if n_20004_0_`n'==1463 & ///
	n_20010_0_`n'<polyp_so_year & n_20010_0_`n'!=-1 & n_20010_0_`n'!=-3
quietly replace polyp_so_age=n_20011_0_`n' if n_20004_0_`n'==1463 & ///
	n_20011_0_`n'<polyp_so_age  & n_20011_0_`n'!=-1 & n_20011_0_`n'!=-3
}

generate polyp_so_time=age_calc-polyp_so_age if polyp_so==1 
replace polyp_so_time=0 if polyp_so_time<0

label variable polyp_so "Self-reported polyp"
label values polyp_so yesno
label variable polyp_so_year "Self-reported polyp procedure: year"
label variable polyp_so_age "Self-reported polyp procedure: age"
label variable polyp_so_time "Self-reported polyp procedure: time since"

*** hospital-reported polyp

generate polyp_h=0
generate polyp_h_date=.

forvalues n=0/212 {
	
quietly replace polyp_h=1 if (s_41270_0_`n'=="K635" | s_41270_0_`n'=="K621") ///
	& ts_41280_0_`n' < basedate
quietly replace polyp_h_date=ts_41280_0_`n' if (s_41270_0_`n'=="K635" | ///
	s_41270_0_`n'=="K621") & ts_41280_0_`n'<basedate & ///
	ts_41280_0_`n'<polyp_h_date
}

generate polyp_h_time=(basedate-polyp_h_date)/365.25 if polyp_h==1 

format polyp_h_date	%td
label variable polyp_h "Hospital-reported polyp"
label values polyp_h yesno
label variable polyp_h_date "Hospital-reported polyp: date"
label variable polyp_h_time "Hospital-reported polyp: time since"

*** hospital-reported polyp procedure
	
generate polyp_ho=0
generate polyp_ho_date=.

forvalues n=0/15 {
	
quietly replace polyp_ho=1 if inlist(substr(s_41272_0_`n',1,3), "H20", ///
	"H23", "H26") & ts_41282_0_`n'<basedate
quietly replace polyp_ho_date=ts_41282_0_`n' if ///
	inlist(substr(s_41272_0_`n',1,3), "H20", "H23", "H26") & ///
	ts_41282_0_`n'<basedate & ts_41282_0_`n'<polyp_ho_date
}

generate polyp_ho_time=(basedate-polyp_ho_date)/365.25 if polyp_ho==1 

format polyp_ho_date %td
label variable polyp_ho "Hospital-reported polyp operation"
label values polyp_ho yesno
label variable polyp_ho_date "Hospital-reported polyp operation: date"
label variable polyp_ho_time "Hospital-reported polyp operation: time since"
	
*** overall polyp before baseline assessment date

generate polyp=0
replace polyp=1 if polyp_s==1 | polyp_so==1 | polyp_h==1  | polyp_ho==1 
label variable polyp "Polyp before baseline"
label values polyp yesno

generate polyp_time=min(polyp_s_time, polyp_so_time, polyp_h_time, ///
	polyp_ho_time) if polyp_s==1 | polyp_so==1 | polyp_h==1 | polyp_ho==1 
label variable polyp_time "Any polyp before baseline: time since"


/* COLORECTAL CANCER SCREENING */

*** self-reported screening (any type)

generate screen_s=n_2345_0_0
recode screen_s (-3 -1 = 0)
label variable screen_s "Screening ever"
label values screen_s yesno

generate screen_s_time=n_2355_0_0
recode screen_s_time (-10 = 0) (-1 -3 = .)
label variable screen_s_time "Screening ever: time since (years)"

*** self-reported screening procedure

generate screen_so=0
generate screen_so_year=.
generate screen_so_age=.

forvalues n=0/31 {
	
quietly replace screen_so=1 if n_20004_0_`n'==1463 | n_20004_0_`n'==1519 
quietly replace screen_so_year=n_20010_0_`n' if (n_20004_0_`n'==1463 | ///
	n_20004_0_`n'==1519) & n_20010_0_`n'<screen_so_year & ///
	n_20010_0_`n'!=-1 & n_20010_0_`n'!=-3
quietly replace screen_so_age=n_20011_0_`n' if (n_20004_0_`n'==1463 | ///
	n_20004_0_`n'==1519) & 	n_20011_0_`n'<screen_so_age & ///
	n_20011_0_`n'!=-1 & n_20011_0_`n'!=-3
}

generate screen_so_time=age_calc-screen_so_age if screen_so==1 
replace screen_so_time=0 if screen_so_time<0

label variable screen_so "Self-reported screening procedure"
label values screen_so yesno
label variable screen_so_year "Self-reported screening procedure: year"
label variable screen_so_age "Self-reported screening procedure: age"
label variable screen_so_time "Self-reported screening procedure: time since"

*** hospital-reported screening procedure before baseline 
	
generate screen_ho=0
generate screen_ho_date=.

forvalues n=0/15 {
	
quietly replace screen_ho=1 if inlist(substr(s_41272_0_`n',1,3), "H20", ///
	"H22", "H23", "H25", "H26", "H28") & ts_41282_0_`n'<basedate
quietly replace screen_ho_date=ts_41282_0_`n' ///
	if inlist(substr(s_41272_0_`n',1,3), "H20", "H22", "H23", "H25", "H26", ///
	"H28") & ts_41282_0_`n'<basedate & ts_41282_0_`n'<screen_ho_date
}

generate screen_ho_time=(basedate-screen_ho_date)/365.25 if screen_ho==1 

format screen_ho_date %td
label variable screen_ho "Hospital-reported screening procedure"
label values screen_ho yesno
label variable screen_ho_date "Hospital-reported screening procedure: date"
label variable screen_ho_time "Hospital-reported screening procedure: time since"

*** screening procedures before baseline assessment date

* any screen
generate screen=0
replace screen=1 if screen_s==1 | screen_so==1 | screen_ho==1 
label variable screen "Any screen before baseline assessment"
label values screen yesno

generate screen_time=min(screen_s_time, screen_so_time, screen_ho_time) if ///
	screen_s==1 | screen_so==1 | screen_ho==1 
label variable screen_time "Any screen before baseline: time"

* screen within last 15 years 

generate screen15=screen
replace screen15=0 if screen_time>15
label variable screen15 "Screening 15 years before baseline"
label values screen15 yesno

generate screen15_time=screen_time
replace screen15_time=. if screen_time>15
label variable screen15_time "Screening 15 years before baseline: time"


* recent screening categories
egen screen15_cat=cut(screen15_time), at(0,2,5,100) icodes
recode screen15_cat (2=1) (1=2) (0=3)
replace screen15_cat=0 if screen15==0
replace screen15_cat=. if screen15==0 & inlist(n_2345_0_0,-1,-3) // no evidence of an answer
label variable screen15_cat "Screening 15 years before baseline: category"
label define screen15_cat_label 0 "Never" 1 "5+ years ago" ///
	2 "2 to <5 years ago" 3 "<2 years ago"
label values screen15_cat screen15_cat_label

* new screening procedure variable (self-reported or from hospital data)
generate screen_new=0
replace screen_new=1 if screen_so==1 | screen_ho==1 
label variable screen_new "Screening procedure before baseline assessment"
label values screen_new yesno

generate screen_new_time=min(screen_so_time, screen_ho_time) if ///
	screen_so==1 | screen_ho==1 
label variable screen_new_time "Screening procedure before baseline: time"
replace screen_new_time=0 if screen_new==0

* new screening procedure variable within 10 years (self-reported or hospital)
generate screen10_new=screen_new
label variable screen10_new "Screening procedure 10 years before baseline assessment"
label values screen10_new yesno
replace screen10_new=0 if screen_new_time>10

generate screen10_new_time=screen_new_time
label variable screen10_new_time "Screening procedure 10 years before baseline: time"
replace screen10_new_time=0 if screen10_new==0


/* INFLAMMATORY BOWEL DISEASE */

* chron's

generate ibd_chrons=0
replace ibd_chrons=1 if ts_131626_0_0<basedate & !inlist(ts_131626_0_0, ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable ibd_chrons "Chron's disease"	
label values ibd_chrons yesno

generate ibd_chrons_time=(basedate-ts_131626_0_0)/365.25 ///
	if ts_131626_0_0<basedate & !inlist(ts_131626_0_0, td(1jan1900), ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable ibd_chrons_time "Time since Chron's disease diagnosed"

* ulcerative colitis

generate ibd_colitis=0
replace ibd_colitis=1 if ts_131628_0_0<basedate & !inlist(ts_131628_0_0, ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable ibd_colitis "Ulcerative colitis"	
label values ibd_colitis yesno	

generate ibd_colitis_time=(basedate-ts_131628_0_0)/365.25 ///
	if ts_131628_0_0<basedate & !inlist(ts_131628_0_0, td(1jan1900), ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable ibd_colitis_time "Time since ulcerative colitis diagnosed"
	
* other IBD	
	
generate ibd_other=0
replace ibd_other=1 if ts_131630_0_0<basedate & !inlist(ts_131630_0_0, ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable ibd_other "Other IBD"	
label values ibd_other yesno	

generate ibd_other_time=(basedate-ts_131630_0_0)/365.25 ///
	if ts_131630_0_0<basedate & !inlist(ts_131630_0_0, td(1jan1900), ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable ibd_other_time "Time since other IBD diagnosed"
	

/* DIABETES */

*** type 1 diabetes

generate diabetes_1=0
replace diabetes_1=1 if ts_130706_0_0<basedate & !inlist(ts_130706_0_0, ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable diabetes_1 "Type 1 diabetes"	
label values diabetes_1 yesno
	
*** type 2 diabetes

generate diabetes_2=0
replace diabetes_2=1 if ts_130708_0_0<basedate & !inlist(ts_130708_0_0, ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable diabetes_2 "Type 2 diabetes"	
label values diabetes_2 yesno

*** unspecified diabetes

generate diabetes_u=0
replace diabetes_u=1 if ts_130714_0_0<basedate & !inlist(ts_130714_0_0, ///
	td(1jan1901), td(2feb1902), td(3mar1903), td(1jan1900))
label variable diabetes_u "Unspecified diabetes"	
label values diabetes_u yesno

*** type 2 or unspecified diabetes_1

generate diabetes=diabetes_2
replace diabetes=1 if diabetes_u==1
label variable diabetes "Type 2 or unspecified diabetes"	
label values diabetes yesno


/* BLOOD MEASUREMENTS */

/*
summ n_30760_0_0 
scalar sd30760=r(sd)
generate n_30760_0_0_sd=n_30760_0_0/sd30760
label variable n_30760_0_0_sd "Per SD of HDL"

summ n_30870_0_0 
scalar sd30870=r(sd)
generate n_30870_0_0_sd=n_30870_0_0/sd30870
label variable n_30870_0_0_sd "Per SD of tryglycerides"

summ n_30780_0_0 
scalar sd30780=r(sd)
generate n_30780_0_0_sd=n_30780_0_0/sd30780
label variable n_30780_0_0_sd "Per SD of LDL"

summ n_30690_0_0 
scalar sd30690=r(sd)
generate n_30690_0_0_sd=n_30690_0_0/sd30690
label variable n_30690_0_0_sd "Per SD of cholesterol"
*/
	
/* MEDICATION */

*** nsaid

generate nsaid=.
replace nsaid=0 if inlist(n_6154_0_0,-7,3,4,5,6)
replace nsaid=1 if n_6154_0_0==1 | n_6154_0_0==2 | n_6154_0_1==2 

label variable nsaid "Taken NSAIDs"
label values nsaid yesno

generate nsaid2=. // this is the response to a question about prescription medication

forvalues n=0/47{

replace nsaid2=1 if inlist(n_20003_0_`n', 1140861806, 1140861808, ///
	1140864860, 1140868226, 1140868282, 1140872040, 1140882108, 1140882190, ///
	1140882192, 1140882268, 1140882392, 1140911760, 1141163138, 1141164044, ///
	1141167844, 1140871310, 1140871374, 1140871388, 1140871394, 1140875540, ///
	1140875616, 1140877962, 1140877964, 1140877966, 1140878030, 1140910496, ///
	1140911086, 1140911748, 1140911750, 1140911762, 1140927152, 1140928656, ///
	1141149110, 1141152166, 1141152168, 1141153082, 1141153134, 1141157412, ///
	1141164254, 1141176278, 1141177836, 1141182814, 1141182868, 1141184226, ///
	1141184546, 1141188652, 1141190952, 1141191742, 1141194296, 1141200576, ///
	1141200748, 1140871462, 1140871472, 1140881612, 1140871168, 1140871174, ///
	1140877892, 1140878034, 1140878036, 1140884488, 1140917394, 1140921828, ///
	1141174424, 1141176878, 1141182674, 1141191028, 1141176662, 1141176668, ///
	1141176670, 1140871542, 1140871546, 1141180140, 1141180148, 1141180150, ///
	1141180152, 1140871336, 1141157452)
}

replace nsaid=1 if nsaid2==1

*** hormone replacement therapy

generate hrt=n_2814_0_0
recode hrt (-1 -3=.)
label variable hrt "HRT ever"
label values hrt yesno

generate hrt_status=hrt
replace hrt_status=2 if hrt==1 & n_3546_0_0==-11
label variable hrt_status "HRT status"
label define hrt_status_label 0 "Never" 1 "Past" 2 "Current"
label values hrt_status hrt_status_label

generate hrt_age_start=n_3536_0_0
recode hrt_age_start (-1 -3=.)
replace hrt_age_start=. if hrt_age_start>age & hrt_age_start!=. // dodgy number
label variable hrt_age_start "Age started HRT (years)"

generate hrt_age_stop=n_3546_0_0 // -11 = still using HRT
recode hrt_age_stop (-1 -3 -11=.)
replace hrt_age_stop=. if n_eid==1101702  // dodgy number
replace hrt_age_stop=44 if n_eid==3045657 // dodgy number
label variable hrt_age_stop "Age stopped HRT (years)"

generate hrt_time=hrt_age_stop-hrt_age_start
replace hrt_time=age-hrt_age_start if n_3546_0_0==-11
label variable hrt_time "Time used hrt (stop-start)"

generate hrt_time_since=age-hrt_age_stop if hrt_status==1
replace hrt_time_since=0 if n_3546_0_0==-11
replace hrt_time_since=. if n_eid==1101702  // dodgy number
label variable hrt_time_since "Time since last used HRT (years)"

egen hrt_time_since_cat=cut(hrt_time_since), at(0,1,5,10,100) icodes
recode hrt_time_since_cat (0=4) (1=3) (3=1) 
replace hrt_time_since_cat=0 if hrt==0 | sex==1
label variable hrt_time_since_cat "Time since last used HRT (category)"
label define hrt_time_since_cat_label 0 "Never" 1 "10 more years" ///
	2 "5 to <10 years" 3 "1 to <5 years" 4 "Current to <1 year"
label values hrt_time_since_cat hrt_time_since_cat_label

*** menopause

generate menop=n_2724_0_0
recode menop (-3=.) (3=2)
replace menop=0 if sex==1
label variable menop "Menopause status"
label define menop_label 0 "No" 1 "Yes" 2 "Unsure"
label values menop menop_label

generate menop_new=menop
label variable menop_new "Menopause status - adjudicated"
* menopausal if unsure and ever took hrt or aged 51 years or older at baseline
replace menop_new=1 if menop==2 & (hrt==1 | age_calc>=51)
replace menop_new=1 if inlist(hrt_time_since_cat,1,2,3,4)
replace menop_new=0 if menop==2 & (hrt==0 | hrt==.) & age_calc<51
* menopausal if ever took HRT and menop=0
replace menop_new=1 if hrt==1 & menop==0
* premenopausal if unsure, never took HRT and <=50 at baseline
replace menop_new=0 if menop==. & hrt==. & sex==0 & age_calc<=50
label values menop_new yesno


*** combined menopausal status and HRT

generate menop_hrt=.
replace menop_hrt=0 if menop_new==0
replace menop_hrt=1 if menop_new==1 & hrt==0
replace menop_hrt=2 if menop_new==1 & hrt==1
label variable menop_hrt "Menopause and HRT"
label define menop_hrt_label 0 "Premenopausal" 1 "Menopause no HRT" ///
	2 "Menopause and HRT"
label values menop_hrt menop_hrt_label


/* EXERCISE */

* code exercise into METs using IPAQ rules

*** temporary mean values used for missing data ***

* walking
generate walk_days=n_864_0_0
recode walk_days (-2=0) (-3 -1=.)
label variable walk_days "Days/week walk >10 mins"

generate walk_dur=n_874_0_0
recode walk_dur (-3 -1=.) (1/9=0) (181/max=180) // truncation rule
label variable walk_dur "Duration of walk >10 mins"

generate walk_pace=n_924_0_0
recode walk_pace (-7 -3=.)
label variable walk_pace "Usual walking pace"
label define walk_pace_label 1 "Slow" 2 "Normal" 3 "Brisk"
label values walk_pace walk_pace_label

generate walk_mets=. // MET minutes per week
replace walk_mets=3.3*walk_days*walk_dur
replace walk_mets=0 if walk_days==0
label variable walk_mets "MET minutes per week - walk"

generate walk_pace_mets=. // MET minutes per week
replace walk_pace_mets=2.5*walk_days*walk_dur if walk_pace==1
replace walk_pace_mets=3.3*walk_days*walk_dur if walk_pace==2 | walk_pace==.
replace walk_pace_mets=5.0*walk_days*walk_dur if walk_pace==3
replace walk_pace_mets=0 if walk_days==0
label variable walk_pace_mets "MET minutes per week - walk (pace adj)"

* moderate activity
generate mod_days=n_884_0_0
recode mod_days (-2=0) (-3 -1=.)
label variable mod_days "Days/week moderate activity >10 mins"

generate mod_dur=n_894_0_0
recode mod_dur (-3 -1=.) (1/9=0) (181/max=180) // truncation rule
label variable mod_dur "Duration of moderate activity >10 mins"

generate mod_mets=. // MET minutes per week
replace mod_mets=4.0*mod_days*mod_dur
replace mod_mets=0 if mod_days==0
label variable mod_mets "MET minutes per week - moderate"

* vigorous activity
generate vig_days=n_904_0_0
recode vig_days (-2=0) (-3 -1=.)
label variable vig_days "Days/week vigorous activity >10 mins"

generate vig_dur=n_914_0_0
recode vig_dur (-3 -1=.) (1/9=0) (181/max=180) // truncation rule
label variable vig_dur "Duration of vigorous activity >10 mins"

generate vig_mets=. // MET minutes per week
replace vig_mets=8.0*vig_days*vig_dur
replace vig_mets=0 if vig_days==0
label variable vig_mets "MET minutes per week - vigorous"

* overall activity
generate activity_mets=(walk_mets + mod_mets + vig_mets)
recode activity_mets (0=1) // for the log transform
label variable activity_mets "MET minutes per week - overall"

generate ln_activity_mets=ln(activity_mets)
label variable ln_activity_mets "ln MET minutes per week"

generate activity_pace_mets=(walk_pace_mets + mod_mets + vig_mets)
recode activity_pace_mets (0=1) // for the log transform
label variable activity_pace_mets "MET minutes per week - overall (pace adj)"

generate ln_activity_pace_mets=ln(activity_pace_mets)
label variable ln_activity_pace_mets "ln MET minutes per week (pace adj)"

*summ ln_activity_mets
*scalar sd_mets=r(sd)
*generate ln_activity_mets_sd=ln_activity_mets/sd_mets
*label variable ln_activity_mets_sd "Per SD of activity METS"


* overall activity category

generate activity_days=walk_days + mod_days + vig_days

generate activity_cat=.
replace activity_cat=0 if activity_mets!=.
replace activity_cat=1 if (inrange(vig_days,3,7) & inrange(vig_dur,20,180)) | ///
	(inrange(mod_days,5,7) & inrange(mod_dur,30,180)) | ///
	(inrange(walk_days,5,7) & inrange(walk_dur,30,180)) 
replace activity_cat=2 if (inrange(vig_days,3,7) & inrange(vig_mets,1500,20000)) | ///
	(inrange(activity_days,7,21) & inrange(activity_mets,3000,20000))
label variable activity_cat "Overall activity category"
label define activity_cat_label 0 "Low" 1 "Moderate" 2 "High"
label values activity_cat activity_cat_label

generate activity_pace_cat=.
replace activity_pace_cat=0 if activity_pace_mets!=.
replace activity_pace_cat=1 if (inrange(vig_days,3,7) & inrange(vig_dur,20,180)) | ///
	(inrange(mod_days,5,7) & inrange(mod_dur,30,180)) | ///
	(inrange(walk_days,5,7) & inrange(walk_dur,30,180)) 
replace activity_pace_cat=2 if (inrange(vig_days,3,7) & inrange(vig_mets,1500,22000)) | ///
	(inrange(activity_days,7,21) & inrange(activity_pace_mets,3000,22000))
label variable activity_pace_cat "Overall activity category (pace adj)"
label values activity_pace_cat activity_cat_label


/* DIET, ALCOHOL AND SMOKING */	

*** calcium supplement

generate supp_calcium=.
replace supp_calcium=0 if inlist(n_6179_0_0,-7,1,2,4,5,6)
replace supp_calcium=1 if inlist(3, n_6179_0_0, n_6179_0_1, n_6179_0_2) 
replace supp_calcium=1 if inlist(7, n_6155_0_0, n_6155_0_1, n_6155_0_2, ///
	n_6155_0_3, n_6155_0_4, n_6155_0_5, n_6155_0_6) // multivitamins
label variable supp_calcium "Take calcium supplement"
label values supp_calcium yesno	

*** fish oil supplement

generate supp_fishoil=.
replace supp_fishoil=0 if inlist(n_6179_0_0,-7,2,3,4,5,6)
replace supp_fishoil=1 if n_6179_0_0==1 
label variable supp_fishoil "Take fish oil supplement"
label values supp_fishoil yesno	

*** vitamin D supplement

generate supp_vitamin_d=.
replace supp_vitamin_d=0 if inlist(n_6155_0_0,-7,1,2,3,5,6)
replace supp_vitamin_d=1 if inlist(4, n_6155_0_0, n_6155_0_1, n_6155_0_2, ///
	n_6155_0_3, n_6155_0_4, n_6155_0_5, n_6155_0_6) | inlist(7, n_6155_0_0, ///
	n_6155_0_1, n_6155_0_2, n_6155_0_3, n_6155_0_4, n_6155_0_5, n_6155_0_6)
label variable supp_vitamin_d "Take vitamin D supplement"
label values supp_vitamin_d yesno	


*** alcohol

generate alcohol=n_1558_0_0
recode alcohol (-3=.) (4 5 6=0) (3=1) (1=3) // 2=2
label variable alcohol "Alcohol use"
label define alcohol_label 0 "Rarely or never" 1 "One or two times per week" ///
	2 "Three or four times per week" 3 "Daily or almost daily"
label values alcohol alcohol_label

generate alcohol_new=alcohol
recode alcohol_new (1=0) (2=1) (3=2)
label define alcohol_new_label 0 "Never or up to two times per week" ///
	1 "Three or four times per week" 2 "Daily or almost daily"
label values alcohol_new alcohol_new_label

*** smoking

generate smoke_ever=n_20116_0_0
recode smoke_ever (-3=.) (2=1)
label variable smoke "Smoking ever"
label values smoke_ever yesno

generate smoke=n_20116_0_0
recode smoke (-3=.)
label variable smoke "Smoking status"
label define smoke_label 0 "Never" 1 "Previous" 2 "Current"
label values smoke smoke_label

generate smoke_time=age-n_2897_0_0 if n_2897_0_0>0
egen smoke_time_cat=cut(smoke_time), at(0,5,20,100) label
replace smoke_time_cat=0 if smoke==2 // current
replace smoke_time_cat=3 if smoke==0 // never
recode smoke_time_cat (3=0) (2=1) (1=2) (0=3)
label variable smoke_time_cat "Smoking time"
label define smoke_time_cat_label 0 "Never" 1 "20+ years ago" ///
	2 "5-19 years ago" 3 "Current or <5 years ago"
label values smoke_time_cat smoke_time_cat_label

*** diet

* oily fish
generate diet_oilyfish=n_1329_0_0 // times per week
recode diet_oilyfish (-10=0) (-1 -3=.) 
label variable diet_oilyfish "Oily fish per week"
label define diet_oilyfish_label 0 "Never" 1 "<1 per week" 2 "1 per week" ///
	3 "2-4 per week" 4 "5-6 per week" 5 "7+ per week"
label values diet_oilyfish diet_oilyfish_label

generate diet_oilyfish_cat=diet_oilyfish
recode diet_oilyfish_cat (4 5=3) 
label variable diet_oilyfish_cat "Oily fish per week: category"
label define diet_oilyfish_cat_label 0 "Never" 1 "<1 per week" ///
	2 "1 per week" 3 "2+ per week"
label values diet_oilyfish_cat diet_oilyfish_cat_label

* oily fish and fish oil supplements
generate fishoil = supp_fishoil
replace fishoil=1 if diet_oilyfish_cat==3
label variable fishoil "Fishoil supplement or oily fish 2+ times per week"
label values fishoil yesno

* processed meat
generate diet_procmeat=n_1349_0_0 // times per week
recode diet_procmeat (-10=0) (-1 -3=.) 
label variable diet_procmeat "Processed meat per week"

generate diet_procmeat_cat=diet_procmeat
recode diet_procmeat_cat (4 5=3) 
label variable diet_procmeat_cat "Processed meat per week: category"
label define diet_procmeat_cat_label 0 "Never" 1 "1 per week" ///
	2 "2 per week" 3 "3+ per week"
label values diet_procmeat_cat diet_procmeat_cat_label

* beef
generate diet_beef=n_1369_0_0 // times per week
recode diet_beef (-10=0) (-1 -3=.) 
label variable diet_beef "Beef per week"

generate diet_beef_cat=diet_beef
recode diet_beef_cat (4 5=3) 
label variable diet_beef_cat "Beef per week: category"
label define diet_beef_cat_label 0 "Never" 1 "1 per week" ///
	2 "2 per week" 3 "3+ per week"
label values diet_beef_cat diet_beef_cat_label

* pork
generate diet_pork=n_1389_0_0 // times per week
recode diet_pork (-10=0) (-1 -3=.) 
label variable diet_pork "Pork per week"

generate diet_pork_cat=diet_pork
recode diet_pork_cat (4 5=3) 
label variable diet_pork_cat "Pork per week: category"
label define diet_pork_cat_label 0 "Never" 1 "1 per week" ///
	2 "2 per week" 3 "3+ per week"
label values diet_pork_cat diet_pork_cat_label

* cereal
generate diet_cereal=n_1458_0_0 // bowls per week
recode diet_cereal (-10=0) (-1 -3=.) (11/max=10) // deal with outliers
label variable diet_cereal "Bowls of cereal per week"

generate diet_cereal_cat=diet_cereal
recode diet_cereal_cat (2 3=1) (4 5 6=2) (7/max=3) 
label variable diet_cereal_cat "Bowls of cereal per week: category"
label define diet_cereal_cat_label 0 "Never" 1 "1-3 per week" ///
	2 "4-6 per week" 3 "7+ per week"
label values diet_cereal_cat diet_cereal_cat_label

generate diet_cereal_cat_new=diet_cereal_cat
recode diet_cereal_cat_new (3=2)
label define diet_cereal_cat_new_label 0 "Never" 1 "1-3 per week" ///
	2 "4+ per week" 
label values diet_cereal_cat_new diet_cereal_cat_new_label
	
generate diet_cereal_fibre_cat=diet_cereal_cat
replace diet_cereal_fibre_cat=0 if n_1468_0_0==5 
label variable diet_cereal_fibre_cat "Bowls of fibre cereal per week: category"
label values diet_cereal_fibre_cat diet_cereal_cat_label

* bread
generate diet_bread=n_1438_0_0 // slices per week
recode diet_bread (-10=0) (-1 -3=.) (51/max=50) // deal with outliers
generate diet_bread_cat=diet_bread
recode diet_bread_cat (2/4=1) (5/10=2) (11/20=3) (21/max=4) 

label define diet_bread_cat_label 0 "None" 1 "1-4 per week" 2 "5-10 per week" ///
	3 "11-20 per week" 4 "21+ per week"
label values diet_bread_cat diet_bread_cat_label

generate diet_bread_whole=diet_bread
replace diet_bread_whole=0 if inlist(n_1448_0_0,1,2,4)  

generate diet_bread_whole_cat=diet_bread_cat
replace diet_bread_whole_cat=0 if inlist(n_1448_0_0,1,2,4)  
label variable diet_bread_whole_cat "Wholemeal/grain bread slices per week: category"
label values diet_bread_whole_cat diet_bread_cat_label

generate diet_bread_white=diet_bread
replace diet_bread_white=0 if inlist(n_1448_0_0,2,3,4)  

generate diet_bread_white_cat=diet_bread_cat
replace diet_bread_white_cat=0 if inlist(n_1448_0_0,2,3,4) 
label variable diet_bread_white_cat "White bread slices per week: category"
label values diet_bread_white_cat diet_bread_cat_label

* cooked vegetables
generate diet_veg_cook=n_1289_0_0
recode diet_veg_cook (-10=0) (-1 -3=.) (11/max=10) // deal with outliers
label variable diet_veg_cook "Cooked vegies per day"

generate diet_veg_cook_cat=diet_veg_cook
recode diet_veg_cook_cat (4/max=3)
label variable diet_veg_cook_cat "Cooked vegies per day: category"
label define diet_veg_cook_cat_label 0 "None" 1 "1 per day" 2 "2 per day" ///
	3 "3+ per day" 
label values diet_veg_cook_cat diet_veg_cook_cat_label

* salad or raw vegetables
generate diet_veg_raw=n_1299_0_0
recode diet_veg_raw (-10=0) (-1 -3=.) (11/max=10) // deal with outliers
label variable diet_veg_raw "Salad or raw vegies per day"

generate diet_veg_raw_cat=diet_veg_raw
recode diet_veg_raw_cat (4/max=3) 
label variable diet_veg_raw_cat "Salad or raw vegies per day: category"
label define diet_veg_raw_cat_label 0 "None" 1 "1 per day" 2 "2 per day" ///
	3 "3+ per day" 
label values diet_veg_raw_cat diet_veg_raw_cat_label

* fresh fruit
generate diet_fruit_fresh=n_1309_0_0
recode diet_fruit_fresh (-10=0) (-1 -3=.) (11/max=10) // deal with outliers
label variable diet_fruit_fresh "Fresh fruit per day"

generate diet_fruit_fresh_cat=diet_fruit_fresh
recode diet_fruit_fresh_cat (4/max=3)
label variable diet_fruit_fresh_cat "Fresh fruit per day: category"
label define diet_fruit_fresh_cat_label 0 "None" 1 "1 per day" 2 "2 per day" ///
	3 "3+ per day" 
label values diet_fruit_fresh_cat diet_fruit_fresh_cat_label

* dried fruit
generate diet_fruit_dried=n_1319_0_0
recode diet_fruit_dried (-10=0) (-1 -3=.)
label variable diet_fruit_dried "Dried fruit per day"

generate diet_fruit_dried_cat=diet_fruit_dried
recode diet_fruit_dried_cat (2/max=1)
label variable diet_fruit_dried_cat "Dried fruit per day: category"
label define diet_fruit_dried_cat_label 0 "None" 1 "1+ per day" 
label values diet_fruit_dried_cat diet_fruit_dried_cat_label


/* CLEAN UP */
	
* drop variables that aren't needed any more (the originally extracted ones)
drop n_20001_0_0-n_20001_0_5 n_20006_0_0-n_20006_0_5 n_20007_0_0-n_20007_0_5 ///
	n_20107_0_0-n_20107_0_9 n_20110_0_0-n_20110_0_10 ///
	n_20111_0_0-n_20111_0_11 s_40002_0_1-s_40002_0_14 ///
	s_40006_0_0-s_40006_16_0 n_40008_0_0-n_40008_16_0 ///
	n_40011_0_0-n_40011_16_0 n_40012_0_0-n_40012_16_0 ///
	s_40013_0_0-s_40013_16_0 ts_40005_0_0-ts_40005_16_0 ///
	n_20002_0_0-n_20002_0_28 n_20004_0_0-n_20004_0_31 ///
	n_20008_0_0-n_20008_0_28 n_20009_0_0-n_20009_0_28 ///
	n_20010_0_0-n_20010_0_31 n_20011_0_0-n_20011_0_31 ///
	s_41270_0_0-s_41270_0_212 s_41271_0_0-s_41271_0_46 ///
	s_41272_0_0-s_41272_0_116 ts_41280_0_0-ts_41280_0_212 ///
	ts_41281_0_0-ts_41281_0_46 ts_41282_0_0-ts_41282_0_116 ///
	n_20003_0_0-n_20003_0_47


