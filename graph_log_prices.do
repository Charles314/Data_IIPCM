*Figure 2

use "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_data_stata.dta", clear

*Get log prices and log shoping trips
gen logprice = log(weighted_price_dep)
gen loga = log(trips)

*Get the average log price paid by shopping trip and good
bysort trips department_code: egen avglogprice=mean(logprice)

*Plot average log price paid across shopping trips
twoway (scatter avglogprice loga if department_code == 1, mcolor(navy) msymbol(t) msize(small)) (scatter avglogprice loga if department_code == 2, mcolor(forest_green) msymbol(d) msize(small)) (scatter avglogprice loga if department_code == 3, mcolor(dkorange) msymbol(o) msize(small)) (scatter avglogprice loga if department_code == 5, mcolor(maroon) msymbol(s) msize(small) legend(ring(0) position(6) size(vsmall) cols(5) label(1 "Dry grocery") label(2 "Frozen foods") label(3 "Dairy") label(4 "Packaged meat") ) xtitle("Log number of shopping trips") ytitle("Average log prices") plotregion(fcolor(white)) graphregion(fcolor(white)))

graph export "C:\Users\charl\OneDrive\Documents\GitHub\graph_price_search.eps", as(eps) name("Graph") preview(off) replace

clear all
