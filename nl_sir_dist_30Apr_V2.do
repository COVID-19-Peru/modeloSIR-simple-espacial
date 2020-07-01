clear all
set more off
set maxvar 32767

cap program drop nlsir
program nlsir, rclass
	version 14.1
	syntax varlist(min=4 max=4) [aw fw iw] if, at(name)
	local inf : word 1 of `varlist'
	local inf_past : word 2 of `varlist'
	local sanos : word 3 of `varlist'
	local densidad : word 4 of `varlist'
	tempname pi beta g delta
	scalar `pi' = `at'[1,1]
	scalar `beta' = `at'[1,2]
	scalar `g' = 0.04
	scalar `delta' = `at'[1,3]
	
	tempvar Wden
	gen double `Wden' = 0 `if'
	forval j = 1($tmax)`=_N' {
		local ubigeo = ubigeo[`j']
		tempvar Wnum`ubigeo'
		gen double `Wnum`ubigeo'' = Poblacion[`j']*exp(-`pi'*DC`ubigeo') `if'
		replace `Wden' = `Wden' + `Wnum`ubigeo'' `if'
	}
	tempvar coef_b
	gen double `coef_b' = 0 `if'
	forval j = 1($tmax)`=_N' {
		local ubigeo = ubigeo[`j']
		tempvar Wi`ubigeo'
		gen double `Wi`ubigeo'' = `Wnum`ubigeo''*`Wden'*i`ubigeo' `if'
		replace `coef_b' = `coef_b' + `Wi`ubigeo'' `if'
	}
	
	replace `inf' = (1-`g')*`inf_past' + `beta'*`coef_b'*`sanos'*(`densidad'^`delta') `if'
end

********************************************************************************
********************************************************************************

** BASE DE DATOS
cap cd D:\Documents\KLV\COVID19
cap cd "\raid\data\KLAB\covid19"
u infectados_MINSA_17Jun, clear
drop if t>55 // Analisis solo hasta el 30 de abril
qui sum t
global tmax = `r(max)' + 1

rename pob Poblacion
drop NuevosInf
tostring ubigeo, replace
replace ubigeo = "0" + ubigeo if strlen(ubigeo)==5
sort ubigeo
by ubigeo: gen InfAcum_past = InfAcum[_n-1]
by ubigeo: gen SanosAcum_past = SanosAcum[_n-1]
by ubigeo: gen int first=1 if _n==1
by ubigeo: gen int last=1 if _n==_N
forval j = 1($tmax)`=_N' {
	local ubigeo = ubigeo[`j']
	qui gen double i`ubigeo' = InfAcum[`j' + t - 1] / Poblacion[`j'] if first!=1
}

* CONSTRUCCION DE DISTANCIAS
** DISTANCIAS EUCLIDEANAS
forval i = 1($tmax)`=_N' {
	local ubigeo = ubigeo[`i']
	g DD`ubigeo' = ((y - y[`i'])^2 + (x - x[`i'])^2)^(1/2)
}
tostring ubigeo, replace
replace ubigeo = "0" + ubigeo if strlen(ubigeo)==5
g ubiprov = substr(ubigeo,1,4)
g ubiprov2 = ubiprov
destring ubiprov, replace
** PEGAR DISTANCIAS PROVINCIALES
preserve
	cap import delimited "C:\Users\Diego\Documents\Covid19\cardistances_prov_corregido.csv", clear encoding("UTF-8")
	cap import delimited "cardistances_prov_corregido.csv", clear encoding("UTF-8")
	quietly {
	cap net install http://www.stata-journal.com/software/sj5-4/dm88_1, replace
	renvars d*, upper
	tempfile distances
	save `distances', replace
	}
restore
merge m:1 ubiprov using `distances', nogen
** DISTANCIAS CORREGIDAS
g DDcapital = .
forval i = 1(1)`=_N' {
	local ubiprov = ubiprov2[`i']
	qui replace DDcapital = DD`ubiprov'01 in `i'
}
forval i = 1($tmax)`=_N' {
	local ubigeo = ubigeo[`i']
	local ubiprov2 = ubiprov2[`i']
	qui g DC`ubigeo' = DDcapital*111 + D`ubiprov2' + DD`ubiprov2'01[`i']*111
}
drop ubiprov2 D0* D1* D2* DD*

* DENSIDAD POBLACIONAL
gen densidad = Poblacion/superficie
gen ldensidad = log(densidad+1)

** ESTIMACION
nl sir @ InfAcum InfAcum_past SanosAcum_past ldensidad if first!=1, parameters(pi beta delta)

** PREDICCION

* IN SAMPLE
quietly {
matrix b = e(b)
local pi = b[1,1]
local beta = b[1,2]
local g = 0.04
local delta = b[1,3]

tempvar Wden
gen double `Wden' = 0 `if'
forval j = 1($tmax)`=_N' {
	local ubigeo = ubigeo[`j']
	tempvar Wnum`ubigeo'
	gen double `Wnum`ubigeo'' = Poblacion[`j']*exp(-`pi'*DC`ubigeo') `if'
	replace `Wden' = `Wden' + `Wnum`ubigeo'' `if'
}
tempvar coef_b
gen double `coef_b' = 0 `if'
forval j = 1($tmax)`=_N' {
	local ubigeo = ubigeo[`j']
	tempvar Wi`ubigeo'
	gen double `Wi`ubigeo'' = `Wnum`ubigeo''*`Wden'*i`ubigeo' `if'
	replace `coef_b' = `coef_b' + `Wi`ubigeo'' `if'
}
}
predictnl infhat = (1-`g')*InfAcum_past + `beta'*`coef_b'*SanosAcum_past*(ldensidad^`delta') if first!=1, se(infhat_se)
drop D* // i0* i1* i2*
save infectados_MINSA_30Apr_predicNL_V2, replace

* OUT OF SAMPLE
*clear all
use infectados_MINSA_30Apr_predicNL_V2, clear
keep if t==55
drop t
rename InfAcum Infectados_t0
rename RecAcum Recuperados_t0
gen Sanos_t0 = Poblacion - Infectados_t0 - Recuperados_t0

qui {
forval t = 1(1)30 {
local t_1 = `t' - 1
	replace `coef_b' = 0
	forval j = 1(1)`=_N' {
		local ubigeo = ubigeo[`j']
		replace i`ubigeo' = Infectados_t`t_1'[`j'] / Poblacion[`j']
		replace `Wi`ubigeo'' = `Wnum`ubigeo''*`Wden'*i`ubigeo'
		replace `coef_b' = `coef_b' + `Wi`ubigeo''
	}
	g Infectados_t`t' = (1-`g')*Infectados_t`t_1' + `beta'*`coef_b'*Sanos_t`t_1'*(ldensidad^`delta')
	g Recuperados_t`t' = Recuperados_t`t_1' + `g'*Infectados_t`t_1'
	g Sanos_t`t' = Poblacion - Infectados_t`t' - Recuperados_t`t'
}
}
cap drop _*
save predicNLInfectados_30Apr30may_V2, replace
