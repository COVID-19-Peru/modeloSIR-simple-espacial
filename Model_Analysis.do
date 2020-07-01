clear all
set more off

global path "D:\Documents\KLV\COVID19"
global path_maps "D:\Documents\Dropbox\Mapas"

u "$path\infectados_MINSA_17Jun", clear
keep ubigeo date t InfAcum RecAcum SanosAcum provincia departamento
rename (InfAcum RecAcum SanosAcum) (O_Infectados O_Recuperados O_Sanos)
g ubigeo2 = ubigeo
tostring ubigeo2, replace
replace ubigeo2 = "0" + ubigeo2 if strlen(ubigeo2)==5
g ubiprov = substr(ubigeo2,1,4)
destring ubiprov, replace
drop ubigeo2
tempfile obs
save `obs', replace

global dates "30Apr30may"
global dates2 "30/04"
global tpers "55"
global Vs "V1 V2 V3"
global top_provs1 "1501 701 2001 1401 1301 2501 2006 401 218"
global top_provs2 "1601 1101 1403 2209 1701 1102 1505 2401 1506"
global top_provs "$top_provs1 $top_provs2"
local i = 1
foreach date in $dates {
	local tper: word `i' of $tpers
	local date2: word `i' of $dates2
	cd "$path\Bases Predic\\`date'"
	foreach V in $Vs {
	if ("`V'"=="V5" | "`V'"=="V6") & ("`date'"=="31Mar30Apr" | "`date'"=="17Jun17Jul") {
		continue
	}
	u predicNLInfectados_`date'_`V', clear
	drop i*
	reshape long Infectados_t Recuperados_t Sanos_t, i(ubigeo) j(t)
	keep ubigeo ubiprov t departamento provincia distrito Infectados_t Recuperados_t Sanos_t
	rename (Infectados_t Recuperados_t Sanos_t) (P_Infectados P_Recuperados P_Sanos)
	replace t = t + `tper'
	destring ubigeo, replace
	merge 1:1 ubigeo t using `obs', nogen
	gen E2 = (O_I - P_I)^2
	sum E2 if t==`tper'+15
	local ECM = round(`r(mean)',3)
	// Nivel nacional
	preserve
		collapse (sum) O_Infectados P_Infectados, by(t)
		label var O_Infectados "Infectados activos (observado)"
		label var P_Infectados "Infectados activos (predicho)"
		label var t "Días tras el primer caso de COVID-19 en Perú"
		drop if t>`tper'+30
		replace P_Infectados=. if P_Infectados<=0
		replace O_Infectados=. if O_Infectados==0
		tw (line O_Infectados t) (line P_Infectados t), xline(`tper')
		graph export "$path\Figuras\\`date'\PlotNac_`date'_`V'.png", as(png) replace
		tw (line O_Infectados t) (line P_Infectados t), yscale(log) xline(`tper')
		graph export "$path\Figuras\\`date'\LPlotNac_`date'_`V'.png", as(png) replace
	restore
	// Nivel top provincial
	preserve
		collapse (sum) O_Infectados P_Infectados (first) provincia departamento, by(ubiprov t)
		label var O_Infectados "Infectados activos (observado)"
		label var P_Infectados "Infectados activos (predicho)"
		label var t "Días tras el primer caso de COVID-19 en Perú"
		label var provincia "Provincia"
		replace provincia=ustrtitle(provincia)
		replace departamento=ustrtitle(departamento)
		gen prov_dep=provincia + ", " + departamento
		label var prov_dep "Provincia, Región"
		drop if t>`tper'+30
		replace P_Infectados=. if P_Infectados<=0
		replace O_Infectados=. if O_Infectados==0
		keep if ubiprov==1501 | ubiprov==701 | ubiprov==2001 | ubiprov==1401 ///
		| ubiprov==1301 | ubiprov==2501 | ubiprov==2006 | ubiprov==401 | ///
		ubiprov==218 | ubiprov==1601 | ubiprov==1101 | ubiprov==1403 | ubiprov==2209 ///
		| ubiprov==1701 | ubiprov==1102 | ubiprov==1505 | ubiprov==2401 | ///
		ubiprov==1506
		tw (line O_Infectados t) (line P_Infectados t), xline(`tper') by(prov_dep, r(6)) ysize(7)
		graph export "$path\Figuras\\`date'\PlotProv_`date'_`V'.png", as(png) replace
		tw (line O_Infectados t) (line P_Infectados t), yscale(log) xline(`tper') by(prov_dep, r(6)) ysize(7)
		graph export "$path\Figuras\\`date'\LPlotProv_`date'_`V'.png", as(png) replace
	restore
	// Nivel bottom provincial
	preserve
		collapse (sum) O_Infectados P_Infectados (first) provincia, by(ubiprov t)
		bys ubiprov (t): g GR = O_Infectados-O_Infectados[_n-15]
		bys ubiprov (t): keep if t==`tper'+15 & O_Infectados[_n-15]<=5
		gsort -GR
		g bottom = 1 in 1/18
		tempfile bottom
		save `bottom', replace
	restore
	preserve
		collapse (sum) O_Infectados P_Infectados (first) provincia departamento, by(ubiprov t)
		label var O_Infectados "Infectados activos (observado)"
		label var P_Infectados "Infectados activos (predicho)"
		label var t "Días tras el primer caso de COVID-19 en Perú"
		label var provincia "Provincia"
		replace provincia=ustrtitle(provincia)
		replace departamento=ustrtitle(departamento)
		gen prov_dep=provincia + ", " + departamento
		label var prov_dep "Provincia, Región"
		drop if t>`tper'+30
		replace P_Infectados=. if P_Infectados<=0
		replace O_Infectados=. if O_Infectados==0
		merge m:1 ubiprov using `bottom', keepus(bottom) nogen
		keep if bottom==1
		tw (line O_Infectados t) (line P_Infectados t), xline(`tper') by(prov_dep, r(6)) ysize(7)
		graph export "$path\Figuras\\`date'\Low_PlotProv_`date'_`V'.png", as(png) replace
		tw (line O_Infectados t) (line P_Infectados t), yscale(log) xline(`tper') by(prov_dep, r(6)) ysize(7)
		graph export "$path\Figuras\\`date'\Low_LPlotProv_`date'_`V'.png", as(png) replace
	restore
	// Mapa Provincial Deciles
	if `i'==5 continue
	preserve
		collapse (sum) O_Infectados P_Infectados (first) departamento provincia, by(ubiprov t)
		label var O_Infectados "Infectados activos (observado)"
		label var P_Infectados "Infectados activos (predicho)"
		label var t "Días tras el primer caso de COVID-19 en Perú"
		label var provincia "Provincia"
		keep if t==`tper'+15
		tostring ubiprov, replace
		replace ubiprov = "0" + ubiprov if strlen(ubiprov)==3
		ren ubiprov IDPROV
		merge 1:1 IDPROV using "$path_maps\data_spmap\prov_d.dta", nogen
		replace departamento=substr(departamento, 1, 1)+ lower(substr(departamento, 2, .))
		replace departamento="Lima y Callao" if regexm(departamento,"Lima") | regexm(departamento,"Callao")
		sum P_Infectados
		local casos = round(r(mean)*196,1)
		spmap P_Infectados using "$path_maps\data_spmap\prov_c.dta", id(_ID) ///
			clmethod(quantile) clnumber(10) fcolor(Heat)					///
			plotregion(margin(medsmall) style(none))						///
			graphregion(margin(zero) style(none))							///
			ocolor(none ..)     ndocolor(none)								///
			title("Predicción") legend(off) saving(P`date'`V', replace)
		spmap O_Infectados using "$path_maps\data_spmap\prov_c.dta", id(_ID) ///
			clmethod(quantile) clnumber(10) fcolor(Heat)					///
			plotregion(margin(medsmall) style(none)) 						///
			graphregion(margin(zero) style(none)) 							///
			ocolor(none ..)     ndocolor(none)								///
			title("Observado") legend(off) saving(O`date'`V', replace)
		graph combine "P`date'`V'" "O`date'`V'",							///
			title("Casos confirmados activos de COVID-19")					///
			note("Mapa coloreado según deciles. Predicción a 15 días con datos al `date2'. ECM=`ECM'.")
		graph export "$path\Figuras\\`date'\MapaProv_`date'_`V'.png", as(png) replace
	restore
	// Mapa Distrital Deciles
	if `i'==5 continue
	preserve
		label var O_Infectados "Infectados activos (observado)"
		label var P_Infectados "Infectados activos (predicho)"
		label var t "Días tras el primer caso de COVID-19 en Perú"
		label var distrito "Distrito"
		keep if t==`tper'+15
		tostring ubigeo, replace
		replace ubigeo = "0" + ubigeo if strlen(ubigeo)==5
		ren ubigeo IDDIST
		merge 1:1 IDDIST using "$path_maps\data_spmap\dist_d.dta", nogen
		replace departamento=substr(departamento, 1, 1)+ lower(substr(departamento, 2, .))
		replace departamento="Lima y Callao" if regexm(departamento,"Lima") | regexm(departamento,"Callao")
		sum P_Infectados
		local casos = round(r(mean)*1874,1)
		spmap P_Infectados using "$path_maps\data_spmap\dist_c.dta", id(_ID) ///
			clmethod(quantile) clnumber(10) fcolor(Heat)					///
			plotregion(margin(medsmall) style(none))						///
			graphregion(margin(zero) style(none))							///
			ocolor(none ..)     ndocolor(none)								///
			title("Predicción") legend(off) saving(P`date'`V', replace)
		spmap O_Infectados using "$path_maps\data_spmap\dist_c.dta", id(_ID) ///
			clmethod(quantile) clnumber(10) fcolor(Heat)					///
			plotregion(margin(medsmall) style(none)) 						///
			graphregion(margin(zero) style(none)) 							///
			ocolor(none ..)     ndocolor(none)								///
			title("Observado") legend(off) saving(O`date'`V', replace)
		graph combine "P`date'`V'" "O`date'`V'",							///
			title("Casos confirmados activos de COVID-19")					///
			note("Mapa coloreado según deciles. Predicción a 15 días con datos al `date2'. ECM=`ECM'.")
		graph export "$path\Figuras\\`date'\MapaDist_`date'_`V'.png", as(png) replace
	restore
	*/
	}
	local ++i
}
