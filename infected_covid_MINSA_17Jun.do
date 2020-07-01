clear all
set more off

cd D:\Documents\KLV\COVID19
import delimited using positivos_covid, clear varn(1)
drop if distrito=="EN INVESTIGACIÓN"
drop if fecha_resultado==""
gen i=1
collapse (sum) infectados=i, by(departamento provincia distrito fecha_resultado)
* MERGE CON UBIGEOS
preserve
	import excel "geodir-ubigeo-inei.xlsx", clear first
	replace Distrito = strupper(Distrito)
	replace Provincia = strupper(Provincia)
	replace Departamento = strupper(Departamento)
	rename (Distrito Provincia Departamento) (distrito provincia departamento)
	replace distrito = subinstr(distrito,"ñ","Ñ",.)
	replace provincia = subinstr(provincia,"ñ","Ñ",.)
	replace departamento = subinstr(departamento,"ñ","Ñ",.)
	replace distrito = subinstr(distrito,"á","Á",.)
	replace provincia = subinstr(provincia,"á","Á",.)
	replace departamento = subinstr(departamento,"á","Á",.)
	replace distrito = subinstr(distrito,"é","É",.)
	replace provincia = subinstr(provincia,"é","É",.)
	replace departamento = subinstr(departamento,"é","É",.)
	replace distrito = subinstr(distrito,"í","Í",.)
	replace provincia = subinstr(provincia,"í","Í",.)
	replace departamento = subinstr(departamento,"í","Í",.)
	replace distrito = subinstr(distrito,"ó","Ó",.)
	replace provincia = subinstr(provincia,"ó","Ó",.)
	replace departamento = subinstr(departamento,"ó","Ó",.)
	replace distrito = subinstr(distrito,"ú","Ú",.)
	replace provincia = subinstr(provincia,"ú","Ú",.)
	replace departamento = subinstr(departamento,"ú","Ú",.)
	replace distrito = subinstr(distrito,"ANDRÉS AVELINO CÁCERES DORREGARAY","ANDRES AVELINO CACERES D.",.)
	replace provincia = subinstr(provincia,"CONTRALMIRANTE VILLA","CONTRALMIRANTE VILLAR",.)
	replace distrito = subinstr(distrito,"CARMEN DE LA LEGUA","CARMEN DE LA LEGUA-REYNOSO",.)
	replace distrito = subinstr(distrito,"26 DE OCTUBRE","VEINTISEIS DE OCTUBRE",.)
	replace distrito = subinstr(distrito,"CORONEL GREGORIO ALBARRACIN LANCHIPA","CORONEL GREGORIO ALBARRACIN L.",.)
	replace distrito = subinstr(distrito,"CONSTITUCIÓN","CONSTITUCION",.)
	replace provincia = subinstr(provincia,"DANIEL ALCIDES CARRI","DANIEL ALCIDES CARRION",.)
	replace distrito = subinstr(distrito,"SAN PEDRO DE PUTINA PUNCO","SAN PEDRO DE PUTINA PUNCU",.)
	replace provincia = subinstr(provincia,"GENERAL SANCHEZ CERR","GENERAL SANCHEZ CERRO",.)
	replace provincia = subinstr(provincia,"ANTONIO RAYMONDI","ANTONIO RAIMONDI",.)
	replace distrito = subinstr(distrito,"JOSÉ MARÍA ARGUEDAS","JOSE MARIA ARGUEDAS",.)
	replace distrito = subinstr(distrito,"RAYMONDI","RAIMONDI",.)
	replace provincia = subinstr(provincia,"CARLOS FERMIN FITZCA","CARLOS FERMIN FITZCARRALD",.)
	replace distrito = subinstr(distrito,"ANCO_HUALLO","ANCO HUALLO",.)
	replace distrito = subinstr(distrito,"DANIEL ALOMIAS ROBLES","DANIEL ALOMIA ROBLES",.)
	replace distrito = subinstr(distrito,"HUAY-HUAY","HUAY HUAY",.)
	replace distrito = subinstr(distrito,"LA YARADA-LOS PALOS","LA YARADA LOS PALOS",.)
	replace distrito = "PAMPAS GRANDE" if distrito=="PAMPAS" & provincia=="HUARAZ"
	replace provincia = "SAN ANTONIO DE PUTINA" if provincia=="SAN ANTONIO DE PUTIN"
	replace distrito = "SAN JUAN DE ISCOS" if distrito=="SAN JUAN DE YSCOS"
	replace provincia = "PUTUMAYO" if provincia=="MAYNAS" & distrito=="PUTUMAYO"
	replace provincia = "PUTUMAYO" if provincia=="MAYNAS" & distrito=="TENIENTE MANUEL CLAVERO"
	replace provincia = "PUTUMAYO" if provincia=="MAYNAS" & distrito=="ROSA PANDURO"
	replace provincia = "PUTUMAYO" if provincia=="MAYNAS" & distrito=="YAGUAS"
	replace distrito = "MARMOT" if distrito=="COMPIN"
	tempfile geodir
	save `geodir', replace
restore
merge m:1 distrito provincia departamento using `geodir', keep(2 3)
rename (Ubigeo Poblacion Superficie Y X) (ubigeo pob superficie y x)
destring ubigeo, replace
replace infectados=0 if infectados==.
drop _merge
g date2 = date(fecha_resultado,"DM19Y") if fecha_resultado!="NULL"
drop fecha_resultado
format date %td
replace date = td(17jun2020) if date==.
tsset ubigeo date
tsfill, full
bysort ubigeo: carryforward departamento provincia distrito pob y x superficie, replace
gsort ubigeo- date
bysort ubigeo: carryforward departamento provincia distrito pob y x superficie, replace
replace infectados = 0 if infectados==.
rename infectados NuevosInf
gen InfAcum = 0
label var InfAcum "Infectados acumulados"
label var NuevosInf "Nuevos infectados"
bysort ubigeo (date): replace InfAcum = sum(NuevosInf)
generate date = string(date2, "%td")
drop date2
bysort ubigeo: gen t = _n - 1
order ubigeo dep prov dist date t NuevosInf InfAcum pob superficie y x
*NUEVOS CASOS
preserve
	drop date InfAcum
	reshape wide NuevosInf, i(ubigeo) j(t)
	order ubigeo dep prov dist pob superficie y x
	save nuevos_inf_covid, replace
	export delimited nuevos_inf_covid.csv, replace
restore
*CASOS ACUMULADOS
preserve
	drop date NuevosInf
	rename InfAcum Infectados_t
	reshape wide Infectados_t, i(ubigeo) j(t)
	order ubigeo dep prov dist pob superficie y x
	save inf_acum_covid, replace
	export delimited inf_acum_covid.csv, replace
restore

* Dias para sanar (SUPUESTO: homogeneo para todos los distritos)
g date2 = date(date,"DM19Y")
format date2 %td
preserve
	import excel PE_Inf_Muertos_Rec_17Jun, first clear
	tempfile heal
	save `heal', replace
restore
merge m:1 date2 using `heal', keepus(days2heal) keep(3)
drop _merge
bys ubigeo (date2) : gen RecAcum = InfAcum[_n-days2heal]
replace RecAcum = 0 if RecAcum==.
bys ubigeo (date2) : gen InfActivos = InfAcum - RecAcum
replace InfActivos = 0 if InfActivos==.
rename InfAcum InfHist
rename InfActivos InfAcum
gen SanosAcum = pob - InfAcum - RecAcum
save infectados_MINSA_17Jun, replace
*CASOS ACUMULADOS (+ RECUPERADOS)
preserve
	drop date date2 days2heal NuevosInf InfHist
	rename (InfAcum RecAcum SanosAcum) (Infectados_t Recuperados_t Sanos_t)
	reshape wide Infectados_t Recuperados_t Sanos_t, i(ubigeo) j(t)
	order ubigeo dep prov dist pob superficie y x
	save inf_acum_covid_Rec, replace
	export delimited inf_acum_covid_Rec.csv, replace
restore
