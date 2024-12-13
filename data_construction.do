*STEP 1: Define the products to be analyzed

import delimited D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\Latest\products.txt

*Focus on nondurables
keep if department_descr == "FROZEN FOODS" | department_descr == "DRY GROCERY" | department_descr == "DAIRY" | department_descr == "PACKAGED MEAT"

*Focus on upcs that are in the Homescan
keep if dataset_found_uc == "ALL" | dataset_found_uc == "HMS"

save "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_products.dta", replace


*STEP 2: Keep purchases in 2011 that are to be analyzed using the products to be analyzed.

import delimited D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\2011\Annual_Files\purchases_2011.csv, clear

*Keep purchases in 2011 that belong to the set of products to be analyzed.
merge m:1 upc upc_ver_uc using D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_products
drop if _merge == 1 | _merge == 2

*Drop a bunch of variables that are not used.
drop _merge upc_desc product_module_descr product_group_descr department_descr brand_descr dataset_found_uc

save "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_purchases_products.dta", replace


*Step 3: Merge trips data with products and compute prices paid

import delimited D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\2011\Annual_Files\trips_2011.csv, clear 

*Join each trip on the information about each of these trips
merge 1:m trip_code_uc using D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_purchases_products

*Drop trips that are not in the set of products to be analyzed.
drop if _merge == 1
drop _merge

*Drop a bunch of useless variables
drop size1_change_flag_uc size1_units size1_amount size1_code_uc retailer_code store_code_uc store_zip3 total_spent deal_flag_uc

*Remove observations (trips) not belonging to the year 2011
generate str1 yearw = ""
replace yearw = substr(purchase_date,1,4)
destring yearw, replace
drop if yearw == 2010 | yearw == 2012
drop yearw

*Create a variable that encodes the month of the trip
gen month = ""
replace month = substr(purchase_date,6,2)
destring month, replace

*Get number of trips for each household month and department
gen dates = 1
egen trips = sum(dates), by(household_code month department_code)

*Get the price paid per unit after coupon discounts
gen price_paid = total_price_paid - coupon_value

*Get the price paid (per unit)
replace price_paid = price_paid/quantity

*Get the weighted price paid within a month
sort household_code month upc upc_ver_uc

egen sum_expenditure = sum(price_paid*quantity), by(household_code month upc upc_ver)
egen sum_quantity = sum(quantity), by(household_code month upc upc_ver)
gen weighted_price = sum_expenditure/sum_quantity

*Average
collapse weighted_price sum_quantity trips, by(household_code month department_code upc upc_ver)

*Get the weighted price paid across 4 department groups
egen sum_expenditure_dep = sum(weighted_price*sum_quantity), by(household_code month department_code)
egen sum_quantity_dep = sum(sum_quantity), by(household_code month department_code)
gen weighted_price_dep = sum_expenditure_dep/sum_quantity_dep

collapse weighted_price_dep sum_quantity_dep trips, by(household_code month department_code)

*Keep months from April to September inclusively
keep if month == 4 | month == 5 | month ==6 | month == 7 | month == 8 | month == 9

save "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_monthly_price_consumption.dta", replace


*Step 4: Merge the data with household information
import delimited D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\2011\Annual_Files\panelists_2011.csv, clear

*Drop a bunch of useless variables
drop fips_state_desc fips_county_desc scantrack_market_identifier_desc dma_name kitchen_appliances tv_items member_1_birth member_1_relationship_sex member_1_employment member_2_birth member_2_relationship_sex member_2_employment member_3_birth member_3_relationship_sex member_3_employment member_4_birth member_4_relationship_sex member_4_employment member_5_birth member_5_relationship_sex member_5_employment member_6_birth member_6_relationship_sex member_6_employment member_7_birth member_7_relationship_sex member_7_employment    

*Need to rename HH_code for linking data sets
rename household_cd household_code

*Put everything together
merge 1:m household_code using D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_monthly_price_consumption
drop if _merge == 1 | _merge == 2
drop _merge

*Drop useless variables
drop projection_factor projection_factor_magnet type_of_residence panelist_zipcd fips_county_cd region_cd scantrack_market_identifier_cd dma_cd household_internet_connection wic_indicator_current wic_indicator_ever_not_current

sort household_code month department_code

*Keep consumers that have a purchase for each good in each month
egen tag = tag( household_code month department_code)
egen total = total(tag), by( household_code )
keep if total == 24

*Keep consumers what are 50 years old and more (don't want online shoppers)
keep if male_head_age >= 7 | female_head_age >= 7

*Keep single households only
keep if household_size == 1

*Remove consumers with negative prices
by household_code: gen negp = sum(tag) if weighted_price_dep <= 0
by household_code: egen totnegp = max(negp)
keep if missing(totnegp)

save "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_data_stata.dta", replace
export delimited sum_quantity_dep weighted_price_dep trips using "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Clean_data\paper_data.csv", replace




******************************************************************************
*Couple Households
******************************************************************************

*Step 4: Merge the data with household information
import delimited D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\2011\Annual_Files\panelists_2011.csv, clear

*Drop a bunch of useless variables
drop fips_state_desc fips_county_desc scantrack_market_identifier_desc dma_name kitchen_appliances tv_items member_1_birth member_1_relationship_sex member_1_employment member_2_birth member_2_relationship_sex member_2_employment member_3_birth member_3_relationship_sex member_3_employment member_4_birth member_4_relationship_sex member_4_employment member_5_birth member_5_relationship_sex member_5_employment member_6_birth member_6_relationship_sex member_6_employment member_7_birth member_7_relationship_sex member_7_employment    

*Need to rename HH_code for linking data sets
rename household_cd household_code

*Put everything together
merge 1:m household_code using D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_monthly_price_consumption
drop if _merge == 1 | _merge == 2
drop _merge

*Drop useless variables
drop projection_factor projection_factor_magnet type_of_residence panelist_zipcd fips_county_cd region_cd scantrack_market_identifier_cd dma_cd household_internet_connection wic_indicator_current wic_indicator_ever_not_current

sort household_code month department_code

*Keep consumers that have a purchase for each good in each month (want people that shop frequently... more susceptible to actually search...)
egen tag = tag( household_code month department_code)
egen total = total(tag), by( household_code )
keep if total == 24

*Keep consumers that are 50 years old and more (don't want online shoppers)
keep if male_head_age >= 7 | female_head_age >= 7

*Keep couple households
keep if household_size == 2

*Remove consumers with negative prices
by household_code: gen negp = sum(tag) if weighted_price_dep <= 0
by household_code: egen totnegp = max(negp)
keep if missing(totnegp)

save "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_data_stata_couples.dta", replace
export delimited sum_quantity_dep weighted_price_dep trips using "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Clean_data\paper_data_couples.csv", replace





******************************************************************************
*Households 3+
******************************************************************************


*Step 4: Merge the data with household information
import delimited D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\2011\Annual_Files\panelists_2011.csv, clear

*Drop a bunch of useless variables
drop fips_state_desc fips_county_desc scantrack_market_identifier_desc dma_name kitchen_appliances tv_items member_1_birth member_1_relationship_sex member_1_employment member_2_birth member_2_relationship_sex member_2_employment member_3_birth member_3_relationship_sex member_3_employment member_4_birth member_4_relationship_sex member_4_employment member_5_birth member_5_relationship_sex member_5_employment member_6_birth member_6_relationship_sex member_6_employment member_7_birth member_7_relationship_sex member_7_employment    

*Need to rename HH_code for linking data sets
rename household_cd household_code

*Put everything together
merge 1:m household_code using D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_monthly_price_consumption
drop if _merge == 1 | _merge == 2
drop _merge

*Drop useless variables
drop projection_factor projection_factor_magnet type_of_residence panelist_zipcd fips_county_cd region_cd scantrack_market_identifier_cd dma_cd household_internet_connection wic_indicator_current wic_indicator_ever_not_current

sort household_code month department_code

*Keep consumers that have a purchase for each good in each month
egen tag = tag( household_code month department_code)
egen total = total(tag), by( household_code )
keep if total == 24

*Keep consumers that are 50 years old and more
keep if male_head_age >= 7 | female_head_age >= 7

*Keep multimember households
keep if household_size > 2

*Remove consumers with negative prices
by household_code: gen negp = sum(tag) if weighted_price_dep <= 0
by household_code: egen totnegp = max(negp)
keep if missing(totnegp)

save "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Master_Files\paper_data_stata_households.dta", replace
export delimited sum_quantity_dep weighted_price_dep trips using "D:\Users\charl\Documents\Nielsen_data\Consumer_panel\Consumer_Panel_Data_2004-2016\Consumer_Panel_Data_2004-2016\nielsen_extracts\HMS\Clean_data\paper_data_households.csv", replace









