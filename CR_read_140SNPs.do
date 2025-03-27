*** UK Biobank -- colorectal cancer dataset extraction -- 140 SNP data ***

version 18.0

cd "C:\Users\gillian.dite\Dropbox\GTG\ColorectalNew\Data\PRS"

import delimited "UKB_Snps_140.txt", varnames(1) clear


rename sample n_eid

* rearrange the SNPs to be in the same order as the Excel file (CRC_140_SNP_list_GSD.xlsx)

order rs72647484 rs4360494 rs12144319 rs7542665 rs6678517 rs17011141 ///
	rs7606562 rs11692435 rs448513 rs11884596 rs983402 rs3731861 rs35470271 ///
	rs9831861 rs6781752 rs13086367 rs12635946 rs72942485 rs10049390 ///
	rs113569514 rs9876206 rs13149359 rs1391441 rs11727676 rs78368589 ///
	rs2735940 rs7708610 rs12514517 rs145364999 rs755229494 rs12659017 ///
	rs4976270 rs2070699 rs1476570 rs3131043 rs116353863 rs116685461 ///
	rs2516420 rs3830041 rs9271695 rs16878812 rs9470361 rs62396735 ///
	rs13204733 rs62404966 rs6928864 rs12672022 rs80077929 rs10951878 ///
	rs3801081 rs16892766 rs6469654 rs117079142 rs6983267 rs7013278 ///
	rs4313119 rs1537372 rs34405347 rs10980628 rs12217641 rs11255841 ///
	rs10821907 rs704017 rs1250567 rs10786560 rs11190164 rs12246635 ///
	rs11196170 rs4450168 rs174533 rs7121958 rs7946853 rs61389091 ///
	rs55864876 rs2186607 rs3087967 rs35808169 rs3217810 rs3217874 ///
	rs10849433 rs2250430 rs2710310 rs77969132 rs11610543 rs12372718 ///
	rs4759277 rs597808 rs1427760 rs7300312 rs9537756 rs7333607 rs45597035 ///
	rs78341008 rs1924816 rs1330889 rs8000189 rs1951864 rs35107139 rs4901473 ///
	rs17094983 rs8020436 rs12708491 rs2293581 rs17816465 rs12594720 ///
	rs56324967 rs745213 rs7495132 rs9924886 rs9930005 rs12447408 rs12149163 ///
	rs62042090 rs4968127 rs73975586 rs1078643 rs983318 rs373585858 ///
	rs75954926 rs11874392 rs34797592 rs28840750 rs1963413 rs12979278 ///
	rs73068325 rs189583 rs994308 rs4813802 rs28488 rs11087784 rs556532366 ///
	rs6058093 rs6031311 rs6066825 rs6067417 rs6063514 rs6091189 rs13831 ///
	rs1741640 rs2738783, after(info)

* labels to use when encoding SNPs

label define A_C_label 0 "C C" 1 "A C" 2 "A A" // A is risk allele
label define C_A_label 0 "A A" 1 "C A" 2 "C C" // C is risk allele

label define A_G_label 0 "G G" 1 "A G" 2 "A A" // A is risk allele
label define G_A_label 0 "A A" 1 "G A" 2 "G G" // G is risk allele

label define A_T_label 0 "T T" 1 "A T" 2 "A A" // A is risk allele
label define T_A_label 0 "A A" 1 "T A" 2 "T T" // T is risk allele

label define C_G_label 0 "G G" 1 "C G" 2 "C C" // C is risk allele
label define G_C_label 0 "C C" 1 "G C" 2 "G G" // G is risk allele

label define C_T_label 0 "T T" 1 "C T" 2 "C C" // C is risk allele
label define T_C_label 0 "C C" 1 "T C" 2 "T T" // T is risk allele

label define G_T_label 0 "T T" 1 "G T" 2 "G G" // G is risk allele
label define T_G_label 0 "G G" 1 "T G" 2 "T T" // T is risk allele

* get rid of missing values (0 0)

foreach v of varlist rs72647484-rs2738783 {

	quietly replace `v'="" if `v'=="0 0"
}

* encode the string genotypes

encode rs72647484, generate(n_rs72647484) label(T_C_label)
encode rs4360494, generate(n_rs4360494) label(G_C_label)
encode rs12144319, generate(n_rs12144319) label(C_T_label)
encode rs7542665, generate(n_rs7542665) label(C_T_label)
encode rs6678517, generate(n_rs6678517) label(A_G_label)
encode rs17011141, generate(n_rs17011141) label(G_A_label)
encode rs7606562, generate(n_rs7606562) label(T_A_label)
encode rs11692435, generate(n_rs11692435) label(G_A_label)
encode rs448513, generate(n_rs448513) label(C_T_label)
encode rs11884596, generate(n_rs11884596) label(C_T_label)
encode rs983402, generate(n_rs983402) label(T_C_label)
encode rs3731861, generate(n_rs3731861) label(T_C_label)
encode rs35470271, generate(n_rs35470271) label(G_A_label)
encode rs9831861, generate(n_rs9831861) label(G_T_label)
encode rs6781752, generate(n_rs6781752) label(A_G_label)
encode rs13086367, generate(n_rs13086367) label(A_G_label)
encode rs12635946, generate(n_rs12635946) label(C_T_label)
encode rs72942485, generate(n_rs72942485) label(G_A_label)
encode rs10049390, generate(n_rs10049390) label(A_G_label)
encode rs113569514, generate(n_rs113569514) label(T_C_label)
encode rs9876206, generate(n_rs9876206) label(C_T_label)
encode rs13149359, generate(n_rs13149359) label(A_C_label)
encode rs1391441, generate(n_rs1391441) label(A_G_label)
encode rs11727676, generate(n_rs11727676) label(C_T_label)
encode rs78368589, generate(n_rs78368589) label(T_C_label)
encode rs2735940, generate(n_rs2735940) label(G_A_label)
encode rs7708610, generate(n_rs7708610) label(A_G_label)
encode rs12514517, generate(n_rs12514517) label(A_G_label)
encode rs145364999, generate(n_rs145364999) label(T_A_label)
encode rs755229494, generate(n_rs755229494) label(G_A_label)
encode rs12659017, generate(n_rs12659017) label(G_A_label)
encode rs4976270, generate(n_rs4976270) label(C_T_label)
encode rs2070699, generate(n_rs2070699) label(T_G_label)
encode rs1476570, generate(n_rs1476570) label(A_G_label)
encode rs3131043, generate(n_rs3131043) label(G_A_label)
encode rs116353863, generate(n_rs116353863) label(C_T_label)
encode rs116685461, generate(n_rs116685461) label(G_A_label)
encode rs2516420, generate(n_rs2516420) label(C_T_label)
encode rs3830041, generate(n_rs3830041) label(T_C_label)
encode rs9271695, generate(n_rs9271695) label(G_A_label)
encode rs16878812, generate(n_rs16878812) label(A_G_label)
encode rs9470361, generate(n_rs9470361) label(A_G_label)
encode rs62396735, generate(n_rs62396735) label(C_T_label)
encode rs13204733, generate(n_rs13204733) label(G_A_label)
encode rs62404966, generate(n_rs62404966) label(C_T_label)
encode rs6928864, generate(n_rs6928864) label(C_A_label)
encode rs12672022, generate(n_rs12672022) label(T_C_label)
encode rs80077929, generate(n_rs80077929) label(T_C_label)
encode rs10951878, generate(n_rs10951878) label(C_T_label)
encode rs3801081, generate(n_rs3801081) label(G_A_label)
encode rs16892766, generate(n_rs16892766) label(C_A_label)
encode rs6469654, generate(n_rs6469654) label(G_C_label)
encode rs117079142, generate(n_rs117079142) label(A_C_label)
encode rs6983267, generate(n_rs6983267) label(G_T_label)
encode rs7013278, generate(n_rs7013278) label(T_C_label)
encode rs4313119, generate(n_rs4313119) label(G_T_label)
encode rs1537372, generate(n_rs1537372) label(G_T_label)
encode rs34405347, generate(n_rs34405347) label(T_G_label)
encode rs10980628, generate(n_rs10980628) label(C_T_label)
encode rs12217641, generate(n_rs12217641) label(C_T_label)
encode rs11255841, generate(n_rs11255841) label(T_A_label)
encode rs10821907, generate(n_rs10821907) label(C_T_label)
encode rs704017, generate(n_rs704017) label(G_A_label)
encode rs1250567, generate(n_rs1250567) label(C_T_label)
encode rs10786560, generate(n_rs10786560) label(G_A_label)
encode rs11190164, generate(n_rs11190164) label(G_A_label)
encode rs12246635, generate(n_rs12246635) label(C_T_label)
encode rs11196170, generate(n_rs11196170) label(A_G_label)
encode rs4450168, generate(n_rs4450168) label(C_A_label)
encode rs174533, generate(n_rs174533) label(G_A_label)
encode rs7121958, generate(n_rs7121958) label(G_T_label)
encode rs7946853, generate(n_rs7946853) label(C_T_label)
encode rs61389091, generate(n_rs61389091) label(C_T_label)
encode rs55864876, generate(n_rs55864876) label(G_A_label)
encode rs2186607, generate(n_rs2186607) label(T_A_label)
encode rs3087967, generate(n_rs3087967) label(T_C_label)
encode rs35808169, generate(n_rs35808169) label(C_T_label)
encode rs3217810, generate(n_rs3217810) label(T_C_label)
encode rs3217874, generate(n_rs3217874) label(T_C_label)
encode rs10849433, generate(n_rs10849433) label(C_T_label)
encode rs2250430, generate(n_rs2250430) label(T_A_label)
encode rs2710310, generate(n_rs2710310) label(C_T_label)
encode rs77969132, generate(n_rs77969132) label(T_C_label)
encode rs11610543, generate(n_rs11610543) label(G_A_label)
encode rs12372718, generate(n_rs12372718) label(G_A_label)
encode rs4759277, generate(n_rs4759277) label(A_C_label)
encode rs597808, generate(n_rs597808) label(G_A_label)
encode rs1427760, generate(n_rs1427760) label(C_T_label)
encode rs7300312, generate(n_rs7300312) label(C_T_label)
encode rs9537756, generate(n_rs9537756) label(C_T_label)
encode rs7333607, generate(n_rs7333607) label(G_A_label)
encode rs45597035, generate(n_rs45597035) label(A_G_label)
encode rs78341008, generate(n_rs78341008) label(C_T_label)
encode rs1924816, generate(n_rs1924816) label(A_G_label)
encode rs1330889, generate(n_rs1330889) label(C_T_label)
encode rs8000189, generate(n_rs8000189) label(T_C_label)
encode rs1951864, generate(n_rs1951864) label(A_G_label)
encode rs35107139, generate(n_rs35107139) label(C_A_label)
encode rs4901473, generate(n_rs4901473) label(G_A_label)
encode rs17094983, generate(n_rs17094983) label(G_A_label)
encode rs8020436, generate(n_rs8020436) label(A_G_label)
encode rs12708491, generate(n_rs12708491) label(G_A_label)
encode rs2293581, generate(n_rs2293581) label(A_G_label)
encode rs17816465, generate(n_rs17816465) label(A_G_label)
encode rs12594720, generate(n_rs12594720) label(C_G_label)
encode rs56324967, generate(n_rs56324967) label(C_T_label)
encode rs745213, generate(n_rs745213) label(G_T_label)
encode rs7495132, generate(n_rs7495132) label(T_C_label)
encode rs9924886, generate(n_rs9924886) label(A_C_label)
encode rs9930005, generate(n_rs9930005) label(C_A_label)
encode rs12447408, generate(n_rs12447408) label(A_G_label)
encode rs12149163, generate(n_rs12149163) label(T_C_label)
encode rs62042090, generate(n_rs62042090) label(T_C_label)
encode rs4968127, generate(n_rs4968127) label(G_A_label)
encode rs73975586, generate(n_rs73975586) label(A_T_label)
encode rs1078643, generate(n_rs1078643) label(A_G_label)
encode rs983318, generate(n_rs983318) label(A_G_label)
encode rs373585858, generate(n_rs373585858) label(A_G_label)
encode rs75954926, generate(n_rs75954926) label(G_A_label)
encode rs11874392, generate(n_rs11874392) label(A_T_label)
encode rs34797592, generate(n_rs34797592) label(T_C_label)
encode rs28840750, generate(n_rs28840750) label(T_G_label)
encode rs1963413, generate(n_rs1963413) label(A_G_label)
encode rs12979278, generate(n_rs12979278) label(T_C_label)
encode rs73068325, generate(n_rs73068325) label(T_C_label)
encode rs189583, generate(n_rs189583) label(G_C_label)
encode rs994308, generate(n_rs994308) label(C_T_label)
encode rs4813802, generate(n_rs4813802) label(G_T_label)
encode rs28488, generate(n_rs28488) label(T_C_label)
encode rs11087784, generate(n_rs11087784) label(G_A_label)
encode rs556532366, generate(n_rs556532366) label(T_C_label)
encode rs6058093, generate(n_rs6058093) label(C_A_label)
encode rs6031311, generate(n_rs6031311) label(T_C_label)
encode rs6066825, generate(n_rs6066825) label(A_G_label)
encode rs6067417, generate(n_rs6067417) label(C_T_label)
encode rs6063514, generate(n_rs6063514) label(C_T_label)
encode rs6091189, generate(n_rs6091189) label(T_C_label)
encode rs13831, generate(n_rs13831) label(G_A_label)
encode rs1741640, generate(n_rs1741640) label(C_T_label)
encode rs2738783, generate(n_rs2738783) label(T_G_label)

* fix the reversed heterozygotes

foreach v of varlist n_rs72647484-n_rs2738783 {

	quietly recode `v' (3=1)
}

* calculate the risk scores for each SNP -- risk=(OR^SNP)/µ

generate r_rs72647484=(1.0517^n_rs72647484)/1.0964
generate r_rs4360494=(1.0386^n_rs4360494)/1.0353
generate r_rs12144319=(1.0683^n_rs12144319)/1.0351
generate r_rs7542665=(1.034^n_rs7542665)/1.0187
generate r_rs6678517=(1.0757^n_rs6678517)/1.0913
generate r_rs17011141=(1.0917^n_rs17011141)/1.0386
generate r_rs7606562=(1.0423^n_rs7606562)/1.07
generate r_rs11692435=(1.0504^n_rs11692435)/1.0928
generate r_rs448513=(1.0054^n_rs448513)/1.0035
generate r_rs11884596=(1.0348^n_rs11884596)/1.0268
generate r_rs983402=(1.0642^n_rs983402)/1.043
generate r_rs3731861=(1.0632^n_rs3731861)/1.0812
generate r_rs35470271=(1.1045^n_rs35470271)/1.0324
generate r_rs9831861=(1.0298^n_rs9831861)/1.0355
generate r_rs6781752=(1.0615^n_rs6781752)/1.0254
generate r_rs13086367=(1.0474^n_rs13086367)/1.0505
generate r_rs12635946=(1.034^n_rs12635946)/1.0426
generate r_rs72942485=(1.056^n_rs72942485)/1.1128
generate r_rs10049390=(1.0466^n_rs10049390)/1.0697
generate r_rs113569514=(1.0423^n_rs113569514)/1.0531
generate r_rs9876206=(1.0463^n_rs9876206)/1.0707
generate r_rs13149359=(1.0534^n_rs13149359)/1.0395
generate r_rs1391441=(1.0149^n_rs1391441)/1.0201
generate r_rs11727676=(1.0093^n_rs11727676)/1.0018
generate r_rs78368589=(1.0818^n_rs78368589)/1.0098
generate r_rs2735940=(1.0904^n_rs2735940)/1.0915
generate r_rs7708610=(1.0391^n_rs7708610)/1.0281
generate r_rs12514517=(1.1066^n_rs12514517)/1.0623
generate r_rs145364999=(1.4185^n_rs145364999)/2.0085
generate r_rs755229494=(1.875^n_rs755229494)/1.0019
generate r_rs12659017=(1.0381^n_rs12659017)/1.0178
generate r_rs4976270=(1.0718^n_rs4976270)/1.0806
generate r_rs2070699=(1.0298^n_rs2070699)/1.0288
generate r_rs1476570=(1.0504^n_rs1476570)/1.0383
generate r_rs3131043=(1.0298^n_rs3131043)/1.0258
generate r_rs116353863=(1.1277^n_rs116353863)/1.0042
generate r_rs116685461=(1.0677^n_rs116685461)/1.1221
generate r_rs2516420=(1.1153^n_rs2516420)/1.225
generate r_rs3830041=(1.0666^n_rs3830041)/1.0187
generate r_rs9271695=(1.093^n_rs9271695)/1.1534
generate r_rs16878812=(1.0809^n_rs16878812)/1.1485
generate r_rs9470361=(1.0555^n_rs9470361)/1.0278
generate r_rs62396735=(1.0336^n_rs62396735)/1.0196
generate r_rs13204733=(1.0664^n_rs13204733)/1.0188
generate r_rs62404966=(1.0751^n_rs62404966)/1.1178
generate r_rs6928864=(1.0545^n_rs6928864)/1.1016
generate r_rs12672022=(1.0067^n_rs12672022)/1.0112
generate r_rs80077929=(1.0093^n_rs80077929)/1.0021
generate r_rs10951878=(1.0545^n_rs10951878)/1.1016
generate r_rs3801081=(1.0256^n_rs3801081)/1.0252
generate r_rs16892766=(1.2336^n_rs16892766)/1.0391
generate r_rs6469654=(1.07^n_rs6469654)/1.0323
generate r_rs117079142=(1.1206^n_rs117079142)/1.0104
generate r_rs6983267=(1.1109^n_rs6983267)/1.1193
generate r_rs7013278=(1.0091^n_rs7013278)/1.0069
generate r_rs4313119=(1.0532^n_rs4313119)/1.0812
generate r_rs1537372=(1.0121^n_rs1537372)/1.0138
generate r_rs34405347=(1.0089^n_rs34405347)/1.0161
generate r_rs10980628=(1.0524^n_rs10980628)/1.0222
generate r_rs12217641=(1.0069^n_rs12217641)/1.0097
generate r_rs11255841=(1.1123^n_rs11255841)/1.1641
generate r_rs10821907=(1.0757^n_rs10821907)/1.1292
generate r_rs704017=(1.0795^n_rs704017)/1.0951
generate r_rs1250567=(1.0481^n_rs1250567)/1.0428
generate r_rs10786560=(1.0082^n_rs10786560)/1.0125
generate r_rs11190164=(1.093^n_rs11190164)/1.0494
generate r_rs12246635=(1.1024^n_rs12246635)/1.0202
generate r_rs11196170=(1.0541^n_rs11196170)/1.0237
generate r_rs4450168=(1.0422^n_rs4450168)/1.0144
generate r_rs174533=(1.0657^n_rs174533)/1.0905
generate r_rs7121958=(1.0811^n_rs7121958)/1.0845
generate r_rs7946853=(1.012^n_rs7946853)/1.0208
generate r_rs61389091=(1.2134^n_rs61389091)/1.452
generate r_rs55864876=(1.0151^n_rs55864876)/1.0279
generate r_rs2186607=(1.0495^n_rs2186607)/1.0519
generate r_rs3087967=(1.1187^n_rs3087967)/1.0703
generate r_rs35808169=(1.0931^n_rs35808169)/1.0323
generate r_rs3217810=(1.1254^n_rs3217810)/1.0317
generate r_rs3217874=(1.0463^n_rs3217874)/1.04
generate r_rs10849433=(1.0479^n_rs10849433)/1.0257
generate r_rs2250430=(1.0615^n_rs2250430)/1.0892
generate r_rs2710310=(1.0146^n_rs2710310)/1.0223
generate r_rs77969132=(1.1715^n_rs77969132)/1.0052
generate r_rs11610543=(1.0485^n_rs11610543)/1.0492
generate r_rs12372718=(1.0937^n_rs12372718)/1.0749
generate r_rs4759277=(1.0289^n_rs4759277)/1.0206
generate r_rs597808=(1.0765^n_rs597808)/1.0806
generate r_rs1427760=(1.0433^n_rs1427760)/1.0461
generate r_rs7300312=(1.0682^n_rs7300312)/1.0795
generate r_rs9537756=(1.0479^n_rs9537756)/1.0595
generate r_rs7333607=(1.0787^n_rs7333607)/1.0373
generate r_rs45597035=(1.0507^n_rs45597035)/1.0671
generate r_rs78341008=(1.011^n_rs78341008)/1.0016
generate r_rs1924816=(1.0519^n_rs1924816)/1.0819
generate r_rs1330889=(1.0463^n_rs1330889)/1.0822
generate r_rs8000189=(1.0484^n_rs8000189)/1.0629
generate r_rs1951864=(1.0059^n_rs1951864)/1.0044
generate r_rs35107139=(1.0955^n_rs35107139)/1.0825
generate r_rs4901473=(1.0476^n_rs4901473)/1.0363
generate r_rs17094983=(1.0062^n_rs17094983)/1.0109
generate r_rs8020436=(1.0298^n_rs8020436)/1.0241
generate r_rs12708491=(1.0475^n_rs12708491)/1.0566
generate r_rs2293581=(1.1329^n_rs2293581)/1.057
generate r_rs17816465=(1.0714^n_rs17816465)/1.0296
generate r_rs12594720=(1.0249^n_rs12594720)/1.0363
generate r_rs56324967=(1.0713^n_rs56324967)/1.0987
generate r_rs745213=(1.0072^n_rs745213)/1.0117
generate r_rs7495132=(1.0463^n_rs7495132)/1.0111
generate r_rs9924886=(1.0565^n_rs9924886)/1.0844
generate r_rs9930005=(1.0061^n_rs9930005)/1.0053
generate r_rs12447408=(1.0079^n_rs12447408)/1.004
generate r_rs12149163=(1.0499^n_rs12149163)/1.0503
generate r_rs62042090=(1.0493^n_rs62042090)/1.0215
generate r_rs4968127=(1.0527^n_rs4968127)/1.0392
generate r_rs73975586=(1.051^n_rs73975586)/1.091
generate r_rs1078643=(1.0776^n_rs1078643)/1.122
generate r_rs983318=(1.0405^n_rs983318)/1.0206
generate r_rs373585858=(1.1166^n_rs373585858)/1.0004
generate r_rs75954926=(1.0922^n_rs75954926)/1.1248
generate r_rs11874392=(1.1742^n_rs11874392)/1.1989
generate r_rs34797592=(1.0859^n_rs34797592)/1.0204
generate r_rs28840750=(1.214^n_rs28840750)/1.4469
generate r_rs1963413=(1.0451^n_rs1963413)/1.056
generate r_rs12979278=(1.0297^n_rs12979278)/1.0317
generate r_rs73068325=(1.0066^n_rs73068325)/1.0024
generate r_rs189583=(1.0827^n_rs189583)/1.0553
generate r_rs994308=(1.0646^n_rs994308)/1.0782
generate r_rs4813802=(1.0853^n_rs4813802)/1.0617
generate r_rs28488=(1.074^n_rs28488)/1.0968
generate r_rs11087784=(1.0913^n_rs11087784)/1.028
generate r_rs556532366=(1.0741^n_rs556532366)/1.0004
generate r_rs6058093=(1.046^n_rs6058093)/1.046
generate r_rs6031311=(1.0369^n_rs6031311)/1.0568
generate r_rs6066825=(1.0745^n_rs6066825)/1.0984
generate r_rs6067417=(1.0337^n_rs6067417)/1.0383
generate r_rs6063514=(1.0562^n_rs6063514)/1.0696
generate r_rs6091189=(1.0564^n_rs6091189)/1.0173
generate r_rs13831=(1.034^n_rs13831)/1.0471
generate r_rs1741640=(1.1214^n_rs1741640)/1.1944
generate r_rs2738783=(1.006^n_rs2738783)/1.0024

* replace missing SNP risks with 1

foreach v of varlist r_rs72647484-r_rs2738783 {

	quietly replace `v'=1 if `v'==.
}

* create PRS (relative risk)

generate prs_rr=1
label variable prs_rr "Colorectal 140-SNP PRS"

foreach v of varlist r_rs72647484-r_rs2738783 {

	quietly replace prs_rr=prs_rr*`v' 
}

generate prs_xb = ln(prs_rr)
label variable prs_xb "Log of colorectal 140-SNP PRS"


* linear combination PRS

generate lc_rs72647484=0.0504*n_rs72647484
generate lc_rs4360494=0.0379*n_rs4360494
generate lc_rs12144319=0.0661*n_rs12144319
generate lc_rs7542665=0.0334*n_rs7542665
generate lc_rs6678517=0.073*n_rs6678517
generate lc_rs17011141=0.0877*n_rs17011141
generate lc_rs7606562=0.0414*n_rs7606562
generate lc_rs11692435=0.0492*n_rs11692435
generate lc_rs448513=0.0054*n_rs448513
generate lc_rs11884596=0.0342*n_rs11884596
generate lc_rs983402=0.0622*n_rs983402
generate lc_rs3731861=0.0613*n_rs3731861
generate lc_rs35470271=0.0994*n_rs35470271
generate lc_rs9831861=0.0294*n_rs9831861
generate lc_rs6781752=0.0597*n_rs6781752
generate lc_rs13086367=0.0463*n_rs13086367
generate lc_rs12635946=0.0334*n_rs12635946
generate lc_rs72942485=0.0545*n_rs72942485
generate lc_rs10049390=0.0455*n_rs10049390
generate lc_rs113569514=0.0414*n_rs113569514
generate lc_rs9876206=0.0453*n_rs9876206
generate lc_rs13149359=0.052*n_rs13149359
generate lc_rs1391441=0.0148*n_rs1391441
generate lc_rs11727676=0.0093*n_rs11727676
generate lc_rs78368589=0.0786*n_rs78368589
generate lc_rs2735940=0.0865*n_rs2735940
generate lc_rs7708610=0.0384*n_rs7708610
generate lc_rs12514517=0.1013*n_rs12514517
generate lc_rs145364999=0.3496*n_rs145364999
generate lc_rs755229494=0.6286*n_rs755229494
generate lc_rs12659017=0.0374*n_rs12659017
generate lc_rs4976270=0.0693*n_rs4976270
generate lc_rs2070699=0.0294*n_rs2070699
generate lc_rs1476570=0.0492*n_rs1476570
generate lc_rs3131043=0.0294*n_rs3131043
generate lc_rs116353863=0.1202*n_rs116353863
generate lc_rs116685461=0.0655*n_rs116685461
generate lc_rs2516420=0.1091*n_rs2516420
generate lc_rs3830041=0.0645*n_rs3830041
generate lc_rs9271695=0.0889*n_rs9271695
generate lc_rs16878812=0.0778*n_rs16878812
generate lc_rs9470361=0.054*n_rs9470361
generate lc_rs62396735=0.033*n_rs62396735
generate lc_rs13204733=0.0643*n_rs13204733
generate lc_rs62404966=0.0724*n_rs62404966
generate lc_rs6928864=0.0531*n_rs6928864
generate lc_rs12672022=0.0067*n_rs12672022
generate lc_rs80077929=0.0093*n_rs80077929
generate lc_rs10951878=0.0531*n_rs10951878
generate lc_rs3801081=0.0253*n_rs3801081
generate lc_rs16892766=0.2099*n_rs16892766
generate lc_rs6469654=0.0677*n_rs6469654
generate lc_rs117079142=0.1139*n_rs117079142
generate lc_rs6983267=0.1052*n_rs6983267
generate lc_rs7013278=0.0091*n_rs7013278
generate lc_rs4313119=0.0518*n_rs4313119
generate lc_rs1537372=0.012*n_rs1537372
generate lc_rs34405347=0.0089*n_rs34405347
generate lc_rs10980628=0.0511*n_rs10980628
generate lc_rs12217641=0.0069*n_rs12217641
generate lc_rs11255841=0.1064*n_rs11255841
generate lc_rs10821907=0.073*n_rs10821907
generate lc_rs704017=0.0765*n_rs704017
generate lc_rs1250567=0.047*n_rs1250567
generate lc_rs10786560=0.0082*n_rs10786560
generate lc_rs11190164=0.0889*n_rs11190164
generate lc_rs12246635=0.0975*n_rs12246635
generate lc_rs11196170=0.0527*n_rs11196170
generate lc_rs4450168=0.0413*n_rs4450168
generate lc_rs174533=0.0636*n_rs174533
generate lc_rs7121958=0.078*n_rs7121958
generate lc_rs7946853=0.0119*n_rs7946853
generate lc_rs61389091=0.1934*n_rs61389091
generate lc_rs55864876=0.015*n_rs55864876
generate lc_rs2186607=0.0483*n_rs2186607
generate lc_rs3087967=0.1122*n_rs3087967
generate lc_rs35808169=0.089*n_rs35808169
generate lc_rs3217810=0.1181*n_rs3217810
generate lc_rs3217874=0.0453*n_rs3217874
generate lc_rs10849433=0.0468*n_rs10849433
generate lc_rs2250430=0.0597*n_rs2250430
generate lc_rs2710310=0.0145*n_rs2710310
generate lc_rs77969132=0.1583*n_rs77969132
generate lc_rs11610543=0.0474*n_rs11610543
generate lc_rs12372718=0.0896*n_rs12372718
generate lc_rs4759277=0.0285*n_rs4759277
generate lc_rs597808=0.0737*n_rs597808
generate lc_rs1427760=0.0424*n_rs1427760
generate lc_rs7300312=0.066*n_rs7300312
generate lc_rs9537756=0.0468*n_rs9537756
generate lc_rs7333607=0.0758*n_rs7333607
generate lc_rs45597035=0.0495*n_rs45597035
generate lc_rs78341008=0.0109*n_rs78341008
generate lc_rs1924816=0.0506*n_rs1924816
generate lc_rs1330889=0.0453*n_rs1330889
generate lc_rs8000189=0.0473*n_rs8000189
generate lc_rs1951864=0.0059*n_rs1951864
generate lc_rs35107139=0.0912*n_rs35107139
generate lc_rs4901473=0.0465*n_rs4901473
generate lc_rs17094983=0.0062*n_rs17094983
generate lc_rs8020436=0.0294*n_rs8020436
generate lc_rs12708491=0.0464*n_rs12708491
generate lc_rs2293581=0.1248*n_rs2293581
generate lc_rs17816465=0.069*n_rs17816465
generate lc_rs12594720=0.0246*n_rs12594720
generate lc_rs56324967=0.0689*n_rs56324967
generate lc_rs745213=0.0072*n_rs745213
generate lc_rs7495132=0.0453*n_rs7495132
generate lc_rs9924886=0.055*n_rs9924886
generate lc_rs9930005=0.0061*n_rs9930005
generate lc_rs12447408=0.0079*n_rs12447408
generate lc_rs12149163=0.0487*n_rs12149163
generate lc_rs62042090=0.0481*n_rs62042090
generate lc_rs4968127=0.0514*n_rs4968127
generate lc_rs73975586=0.0497*n_rs73975586
generate lc_rs1078643=0.0747*n_rs1078643
generate lc_rs983318=0.0397*n_rs983318
generate lc_rs373585858=0.1103*n_rs373585858
generate lc_rs75954926=0.0882*n_rs75954926
generate lc_rs11874392=0.1606*n_rs11874392
generate lc_rs34797592=0.0824*n_rs34797592
generate lc_rs28840750=0.1939*n_rs28840750
generate lc_rs1963413=0.0441*n_rs1963413
generate lc_rs12979278=0.0293*n_rs12979278
generate lc_rs73068325=0.0066*n_rs73068325
generate lc_rs189583=0.0795*n_rs189583
generate lc_rs994308=0.0626*n_rs994308
generate lc_rs4813802=0.0819*n_rs4813802
generate lc_rs28488=0.0714*n_rs28488
generate lc_rs11087784=0.0874*n_rs11087784
generate lc_rs556532366=0.0715*n_rs556532366
generate lc_rs6058093=0.045*n_rs6058093
generate lc_rs6031311=0.0362*n_rs6031311
generate lc_rs6066825=0.0719*n_rs6066825
generate lc_rs6067417=0.0331*n_rs6067417
generate lc_rs6063514=0.0547*n_rs6063514
generate lc_rs6091189=0.0549*n_rs6091189
generate lc_rs13831=0.0334*n_rs13831
generate lc_rs1741640=0.1146*n_rs1741640
generate lc_rs2738783=0.006*n_rs2738783

* replace missing SNP risks with 1

foreach v of varlist lc_rs72647484-lc_rs2738783 {

	quietly replace `v'=0 if `v'==.
}

* create PRS (linear combination)

generate lc_prs_xb=0
label variable lc_prs_xb "Colorectal 140-SNP PRS linear combination"

foreach v of varlist lc_rs72647484-lc_rs2738783 {

	quietly replace lc_prs_xb=lc_prs_xb+`v' 
}



* count genotyped SNPs

generate snp_n=0
label variable snp_n "Number of SNPs genotyped"

foreach v of varlist n_rs72647484-n_rs2738783 {

	quietly replace snp_n=snp_n+1 if `v'!=.
}

save "CR_140PRS.dta", replace

cd ..

keep n_eid sex prs_rr prs_xb lc_prs_xb snp_n 

save "CR_140PRS_risks.dta", replace
