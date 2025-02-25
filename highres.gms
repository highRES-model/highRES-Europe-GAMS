******************************************
* highres main script
******************************************

* $ontext
option profile=1
* $offtext
option limrow=0, limcol=0, solprint=OFF
option decimals = 4
$offlisting
$ONMULTI
$ONEPS
$offdigit

* Switches:

* log = text file to store details about model run time, optimality or not, etc.
* gdx2sql (ON/OFF) = whether to convert output GDX to sqlite database -> easier
*                       to read into Python
* outname = output name of GDX file
*
* GWatts (YES/NO) = model is run in GW (YES) or MW (NO)
*
* storage (ON/OFF) = should storage be included in the model run
* hydrores (ON/OFF) = should reservoir hydro be incliuded in the model run
* sensitivity (ON/OFF) = whether a sensitivity file is available
* UC (ON/OFF) = unit committment switch
* store_uc (ON/OFF) = storage unit committment switch
*
* f_res (ON/OFF) = should frequency response requirements be modelled
*
* sense_run = sensitivity file identifier
* esys_scen = energy system scenario (sets the carbon budget and demands to be
*               used)
* psys_scen = power system scenario (sets which technologies are available)
* RPS = renewable portfolio standard
* vre_restrict = VRE land use deployment scenario name
* model_yr = which year in the future are we modelling
* weather_yr = which weather year do we use
* dem_yr = which demand year do we use
* fx_trans (YES/NO) = fix transmission network to input values
* fx_caps_to = file containing capacities to fix the system to
* co2intensity = Average annual carbon dioxide emission intensity
*                the model is allowed to produce (gCO2/kWh)
* emis_price = emission price for non VRE
* store_initial_level = energy stored at the start (%/100)
* store_final_level = energy stored at the end (%/100)
* hydro_res_initial_fill = how full are the reservoirs at the start and end (%/100)
* hydro_res_min = minimum reservoir level (%/100)
*
* pen_gen (ON/OFF) = whether value of lost load (VoLL) is modelled
* pgen = pen_gen cost
*
*
* DISABLED switches
*
* fx_natcap (YES/NO) = fix total national capacities ->  let highRES decide
*                          where to place them
* water (ON/OFF) = model technologies with a water footprint
*
*
* Paths, logging and IO

$setglobal datafolderpath "."
$setglobal codefolderpath "."
$setglobal log "test_log"
$setglobal gdx2sql "OFF"
$setglobal outname "results"

* Units

$setglobal GWatts "YES"

* Modules

$setglobal storage "ON"
$setglobal hydrores "ON"
$setglobal sensitivity "OFF"
$setglobal UC "OFF"
$setglobal store_uc "OFF"

** Unit comittment switches

$setglobal f_res "OFF"

* Case options

$setglobal sense_run "None"
$setglobal esys_scen "BASE"
$setglobal psys_scen "BASE"
$setglobal RPS "optimal"
$setglobal vre_restrict ""
$setglobal model_yr "2050"
$setglobal weather_yr "2010"
$setglobal dem_yr "2010"
$setglobal fx_trans "NO"
$setglobal fx_caps_to ""
$setglobal co2intensity "2"
$setglobal emis_price "0"

$setglobal store_initial_level "0.5"
$setglobal store_final_level "0.5"

$setglobal hydro_res_initial_fill "0.8"
$setglobal hydro_res_min "0.5"

$set pen_gen "OFF"
$setglobal pgen "20.0"

* Disabled switches
$setglobal water "OFF"
$setglobal fx_natcap "NO"



**************************************************


* rescale from MW to GW for better numerics (allegedly)

scalar MWtoGW;

$ifThen "%GWatts%" == YES

MWtoGW=1E3;

$else

MWtoGW=1;

$endif

$INCLUDE %codefolderpath%/highres_data_input.gms


$IF "%storage%" == ON $INCLUDE %codefolderpath%/highres_storage_setup.gms

* WARNING: for parameter updates to work there can be no arithmetic in the code
* before the update is run -> sensitivity data must be imported here


$IF "%sensitivity%" == ON $INCLUDE %datafolderpath%/sensitivity_%sense_run%.dd

* if no RPS set just do an optimal run

$IF "%RPS%" == "optimal" $GOTO optimal1

scalar
RPS
/%RPS%/
;
RPS=RPS/100.

$label optimal1


demand(z,h)=demand(z,h)/MWtoGW;
gen_cap2area(vre)=gen_cap2area(vre)/MWtoGW;
trans_links_cap(z,z_alias,trans)=trans_links_cap(z,z_alias,trans)/MWtoGW;
gen_unitsize(non_vre)=gen_unitsize(non_vre)/MWtoGW;
gen_maxramp(non_vre)=gen_maxramp(non_vre)/MWtoGW;




* Existing VRE capacity aggregated to zones

* exist_vre_cap_r(vre,z,r) = 0.0;

$ontext
gen_exist_pcap_z(z,vre,"FX")=sum(r,exist_vre_cap_r(vre,z,r));

$offtext

* Existing zonal capacity aggregated to national

parameter gen_exist_cap(g);
gen_exist_cap(g)=sum((z,lt),gen_exist_pcap_z(z,g,lt));

* Limit which regions a given VRE tech can be built in
* based on buildable area in that region. Stops offshore solar
* and onshore offshore wind.

set vre_lim(vre,z,r);
vre_lim(vre,z,r)=((area(vre,z,r)+sum(lt,gen_exist_pcap_r(vre,z,r,lt)))>0.);

* Non VRE cap lim to dynamic set, stops Nuclear being built in certain countries
*   (e.g. Austria)

set gen_lim(z,g);
gen_lim(z,non_vre)=((sum(lt,gen_lim_pcap_z(z,non_vre,lt))
    +sum(lt,gen_exist_pcap_z(z,non_vre,lt)))>0.);
gen_lim(z,vre)=(sum(r,(area(vre,z,r)+sum(lt,gen_exist_pcap_r(vre,z,r,lt))))>0.);


sets
gen_lin(z,non_vre)
ramp_on(z,non_vre)
mingen_on(z,non_vre)
ramp_and_mingen(z,non_vre)
;



$ifThen "%UC%" == ON



* set for generations that can provide quick start operating reserve

set gen_quick(non_vre);

* only OCGT can provide quick start reserves

gen_quick("NaturalgasOCGTnew")=YES;

* generators that are represented as continous linear capacity chunks

gen_lin(z,non_vre)=not ((gen_uc_lin(non_vre) or gen_uc_int(non_vre)) and
    uc_z(z));

* generators that are continous linear chunks can have a mingen of 0

parameter gen_mingen_lin(non_vre);
gen_mingen_lin(non_vre)=0.;

$else

* if UC is not on all generators are linear chunks

gen_lin(z,non_vre)=YES;

parameter gen_mingen_lin(non_vre);
gen_mingen_lin(non_vre)=gen_mingen(non_vre);
gen_mingen_lin("NaturalgasOCGTnew")=0.0;
gen_mingen_lin("NaturalgasCCGTwithCCSnewOT")=0.0;

$endIf

* Sets to ensure ramp/mingen constraints are only created where relevant

ramp_on(z,non_vre)=((gen_maxramp(non_vre)*60./gen_unitsize(non_vre)) < 1.0 and
    gen_lim(z,non_vre) and gen_lin(z,non_vre) and gen_unitsize(non_vre) > 0.);

mingen_on(z,non_vre)=(gen_mingen_lin(non_vre) > 0. and gen_lim(z,non_vre) and
    gen_lin(z,non_vre));

ramp_and_mingen(z,non_vre) = (ramp_on(z,non_vre) or mingen_on(z,non_vre));

* Buildable area per cell from km2 to MW power capacity

area(vre,z,r)=area(vre,z,r)$(vre_lim(vre,z,r))*gen_cap2area(vre);

* To be conservative, existing capacity is removed from new capacity limit

area(vre,z,r)=area(vre,z,r)-sum(lt,gen_exist_pcap_r(vre,z,r,lt));
area(vre,z,r)$(area(vre,z,r)<0.) = 0.;

* Fuel, varom and emission costs for non VRE gens;

gen_varom(non_vre)=gen_fuelcost(non_vre)+gen_emisfac(non_vre)*%emis_price%+
    gen_varom(non_vre);


* Penalty generation setup
* VoLL set at 20000£/MWh


* Solar marginal - small value necessary to avoid transmission system issue

gen_varom("Solar")=0.001;

*trans_varom(trans)=0.001;

* Rescale parameters for runs that are greater or less than one year

if (card(h) < 8760,
*co2_budget=round(co2_budget*(card(h)/8760.),8);
gen_capex(g)=round(gen_capex(g)*(card(h)/8760.),8);
gen_fom(g)=round(gen_fom(g)*(card(h)/8760.),8);
trans_line_capex(trans)=round(trans_line_capex(trans)*(card(h)/8760.),8);
trans_sub_capex(trans)=round(trans_sub_capex(trans)*(card(h)/8760.),8);

store_fom(s)=round(store_fom(s)*(card(h)/8760.),8);
store_pcapex(s)=round(store_pcapex(s)*(card(h)/8760.),8);
store_ecapex(s)=round(store_ecapex(s)*(card(h)/8760.),8);
);




Variables
costs                                    total electricty system costs

* Total cost components

costs_gen_capex(z)
costs_gen_fom(z)
costs_gen_varom(z)
$IF "%UC%" == ON costs_gen_start(z)
costs_store_capex(z)
costs_store_fom(z)
costs_store_varom(z)
$IF "%store_uc%" == ON costs_store_start(z)
costs_trans_capex(z)
costs_trans_fom(z)
$IF "%pen_gen%" == ON costs_pgen(z)

Positive variables
var_new_pcap(g)                    new generation capacity at national level
var_new_pcap_z(z,g)                new generation capacity at zonal level
var_exist_pcap(g)                  existing generation capacity at national
*                                  level
var_exist_pcap_z(z,g)              existing generation capacity at zonal level
var_tot_pcap(g)                    total generation capacity at national level
var_tot_pcap_z(z,g)                total generation capacity at zonal level
var_gen(h,z,g)                     generation by hour and technology
var_new_vre_pcap_r(z,vre,r)        new VRE capacity at grid cell level by
*                                  technology and zone
var_exist_vre_pcap_r(z,vre,r)      existing VRE capacity at grid cell level by
*                                  technology and zone
var_vre_gen_r(h,z,vre,r)           VRE generation at grid cell level by hour
*                                  zone and technology
var_vre_curtail(h,z,vre,r)         VRE power curtailed
*var_non_vre_curtail(z,h,non_vre)
var_trans_flow(h,z,z_alias,trans)  Flow of electricity from node to node by hour
*                                  (MW)
var_trans_pcap(z,z_alias,trans)    Capacity of node to node transmission links
*                                  (MW)

var_pgen(h,z)                      Penalty generation

;

* Synchronous condensers can have negative generation - they require energy to
*   function

* var_gen.LO(h,z,"SynCon")=-inf;

*** Transmission set up ***

* Sets up bidirectionality of links

trans_links(z_alias,z,trans)$(trans_links(z,z_alias,trans))
    =trans_links(z,z_alias,trans);

trans_links_cap(z_alias,z,trans)$(trans_links_cap(z,z_alias,trans) > 0.)
    =trans_links_cap(z,z_alias,trans);

trans_links_dist(z,z_alias,trans)=trans_links_dist(z,z_alias,trans)/100.;

* Bidirectionality of link distances for import flow reduction -> both monodir
*   and bidir needed, former for capex.

parameter trans_links_dist_bidir(z,z_alias,trans);

trans_links_dist_bidir(z,z_alias,trans)=trans_links_dist(z,z_alias,trans);
trans_links_dist_bidir(z_alias,z,trans)$(trans_links_dist(z,z_alias,trans) > 0.)
    =trans_links_dist(z,z_alias,trans);

* Set transmission capacities to historic

$ifThen "%fx_trans%" == "YES"

var_trans_pcap.FX(z,z_alias,trans)$(trans_links(z,z_alias,trans))
    = trans_links_cap(z,z_alias,trans);

$else

* Or limit all links to some maximum

var_trans_pcap.UP(z,z_alias,trans)$(trans_links(z,z_alias,trans))=50.;

$endIf



*******************************

var_exist_pcap_z.UP(z,g)$(gen_exist_pcap_z(z,g,"UP"))
    = gen_exist_pcap_z(z,g,"UP");
*var_exist_pcap_z.L(z,g)$(gen_exist_pcap_z(z,g,"UP"))
*    = gen_exist_pcap_z(z,g,"UP");

var_exist_pcap_z.LO(z,g)$(gen_exist_pcap_z(z,g,"LO"))
    = gen_exist_pcap_z(z,g,"LO");
var_exist_pcap_z.FX(z,g)$(gen_exist_pcap_z(z,g,"FX"))
    = gen_exist_pcap_z(z,g,"FX");

var_exist_pcap_z.UP(z,g)$(not (sum(lt,gen_exist_pcap_z(z,g,lt)) > 0.)) = 0.0;



var_tot_pcap_z.UP(z,g)$(gen_lim_pcap_z(z,g,'UP'))=gen_lim_pcap_z(z,g,'UP');
var_tot_pcap_z.LO(z,g)$(gen_lim_pcap_z(z,g,'LO'))=gen_lim_pcap_z(z,g,'LO');
var_tot_pcap_z.FX(z,g)$(gen_lim_pcap_z(z,g,'FX'))=gen_lim_pcap_z(z,g,'FX');


var_exist_vre_pcap_r.UP(z,vre,r)$(gen_exist_pcap_r(vre,z,r,"UP"))=gen_exist_pcap_r(vre,z,r,"UP");
var_exist_vre_pcap_r.LO(z,vre,r)$(gen_exist_pcap_r(vre,z,r,"LO"))=gen_exist_pcap_r(vre,z,r,"LO");
var_exist_vre_pcap_r.FX(z,vre,r)$(gen_exist_pcap_r(vre,z,r,"FX"))=gen_exist_pcap_r(vre,z,r,"FX");

* turn off the potential for any existing capacity that is not given to the model

parameter gen_exist_r_on(z,g);
gen_exist_r_on(z,vre)=sum((r,lt),gen_exist_pcap_r(vre,z,r,lt));

var_exist_pcap_z.UP(z,g)$(not (sum(lt,gen_exist_pcap_z(z,g,lt)) > 0. or gen_exist_r_on(z,g) > 0.)) = 0.0;

* $IF "%fx_natcap%" == YES var_new_pcap.FX(g)$(gen_fx_natcap(g))=gen_fx_natcap(g);



$ifThen "%UC%" == ON
* parameters for UC, used by both gen and storage UC so needs to be defined
*   here

    scalars
    f_res_time      frequency response ramp up window (minutes)    / 0.083 /
    res_time        operating reserve ramp up window (minutes)     / 20. /
    unit_cap_lim_z  limit the maximum capacity of each units deployed in each
        zone (MW)                                                  /50000./
    res_margin      operating reserve margin (fraction of demand)  / 0.1 /
    ;

    sets
    service_type    ancillary service type                / f_response, reserve/
    hh_minup(h)     max minup hours                       / 0*23 /
    hh_mindown(h)                                         / 0*7 /;
    ;
$endIf

$IF "%store_uc%" == ON $INCLUDE %codefolderpath%/highres_storage_uc_setup.gms

$IF "%UC%" == ON $INCLUDE %codefolderpath%/highres_uc_setup.gms

*$IF "%water%" == ON $INCLUDE %codefolderpath%/highres_water_setup.gms

$IF "%hydrores%" == ON $INCLUDE %codefolderpath%/highres_hydro.gms



$ifThen NOT "%fx_caps_to%" == ""

parameters
par_new_pcap_z(z,g)
par_exist_pcap_z(z,g)
par_new_store_pcap_z(z,s)
par_exist_store_pcap_z(z,s)
par_new_store_ecap_z(z,s)
par_exist_store_ecap_z(z,s)
par_trans_pcap(z,z_alias,trans)
;

$INCLUDE %datafolderpath%/%fx_caps_to%.dd
;

var_new_pcap_z.FX(z,g) = par_new_pcap_z(z,g);
var_exist_pcap_z.FX(z,g)=par_exist_pcap_z(z,g);

* redefine area, vre_lim and gen_lim so all VREs can be fixed

parameter area(vre,z,r);
area(vre,z,r)$(ord(z) eq ord(r))=par_new_pcap_z(z,vre);

vre_lim(vre,z,r)=((area(vre,z,r)+exist_vre_cap_r(vre,z,r))>0.);
gen_lim(z,vre)=sum(r,(area(vre,z,r)>0.));

var_new_store_pcap_z.FX(z,s)=par_new_store_pcap_z(z,s);
var_exist_store_pcap_z.FX(z,s)=par_exist_store_pcap_z(z,s);
var_new_store_ecap_z.FX(z,s)=par_new_store_ecap_z(z,s);
var_exist_store_ecap_z.FX(z,s)=par_exist_store_ecap_z(z,s);

var_trans_pcap.FX(z,z_alias,trans)=par_trans_pcap(z,z_alias,trans);


$endIf


Equations
eq_obj

eq_costs_gen_capex
eq_costs_gen_fom
eq_costs_gen_varom
$IF "%UC%" == ON eq_costs_gen_start
eq_costs_store_capex
eq_costs_store_fom
eq_costs_store_varom
$IF "%store_uc%" == ON eq_costs_store_start
eq_costs_trans_capex
eq_costs_trans_fom

eq_elc_balance

eq_new_pcap
eq_exist_pcap
eq_tot_pcap
eq_tot_pcap_z

eq_gen_max
eq_gen_min
eq_ramp_up
eq_ramp_down
*eq_curtail_max_non_vre

eq_new_vre_pcap_z
eq_exist_vre_pcap_z
eq_gen_vre
eq_gen_vre_r

eq_area_max

eq_trans_flow
eq_trans_bidirect
$IF "%pen_gen%" == ON eq_pen_gen

eq_co2_target

*eq_cap_margin

;

******************************************
* OBJECTIVE FUNCTION
******************************************

eq_obj .. costs =E= sum(z,

costs_gen_capex(z)
+costs_gen_fom(z)
+costs_gen_varom(z)
$IF "%UC%" == ON +costs_gen_start(z)
+costs_store_capex(z)
+costs_store_fom(z)
+costs_store_varom(z)
$IF "%store_uc%" == ON +costs_store_start(z)
$IF "%pen_gen%" == ON +costs_pgen(z)
+costs_trans_capex(z)
+costs_trans_fom(z));


eq_costs_gen_capex(z)..
    costs_gen_capex(z) =E= sum(g,var_new_pcap_z(z,g)*gen_capex(g));

eq_costs_gen_fom(z)..
    costs_gen_fom(z) =E= sum(g,var_new_pcap_z(z,g)*gen_fom(g))
                            +sum(g,var_exist_pcap_z(z,g)*gen_fom(g));

eq_costs_gen_varom(z)..
    costs_gen_varom(z) =E= sum((h,gen_lim(z,g)),var_gen(h,z,g)*gen_varom(g));

$ifThen "%UC%" == ON
    eq_costs_gen_start(z)..
        costs_gen_start(z) =E= sum((h,non_vre)$(gen_uc_int(non_vre) and
            gen_lim(z,non_vre)),var_up_units(h,z,non_vre)
            *gen_startupcost(non_vre))
        +sum((h,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)),
            var_up_units_lin(h,z,non_vre)*gen_startupcost(non_vre));
$endIf

$ifThen "%pen_gen%" == ON
    eq_pen_gen(z)..
        costs_pgen(z) =E= sum((h),var_pgen(h,z)*%pgen%);
$endIf

eq_costs_store_capex(z)..
    costs_store_capex(z) =E= sum(s,var_new_store_pcap_z(z,s)*store_pcapex(s)
        +var_new_store_ecap_z(z,s)*store_ecapex(s));

eq_costs_store_fom(z)..
    costs_store_fom(z) =E= sum(s,var_exist_store_pcap_z(z,s)*store_fom(s))
        +sum(s,var_new_store_pcap_z(z,s)*store_fom(s));

eq_costs_store_varom(z)..
    costs_store_varom(z) =E= sum((h,s)$s_lim(z,s),var_store_gen(h,z,s)
        *store_varom(s));

$IF "%store_uc%" == ON eq_costs_store_start(z)..
$IF "%store_uc%" == ON     costs_store_start(z) =E= sum(
$IF "%store_uc%" == ON          (h,s_lim(z,store_uc_lin)),
$IF "%store_uc%" == ON          var_store_up_units_lin(h,z,store_uc_lin)
$IF "%store_uc%" == ON          *store_startupcost(store_uc_lin));

eq_costs_trans_capex(z)..
    costs_trans_capex(z) =E= sum(trans_links(z,z_alias,trans),
        var_trans_pcap(z,z_alias,trans)*trans_links_dist(z,z_alias,trans)
        *trans_line_capex(trans))
    +sum(trans_links(z,z_alias,trans),
        var_trans_pcap(z,z_alias,trans)$(trans_links_dist(z,z_alias,trans))
        *trans_sub_capex(trans)*2);

* assume 2% fom costs for transmission

eq_costs_trans_fom(z) .. costs_trans_fom(z) =E= costs_trans_capex(z)*0.02;



******************************************
* SUPPLY-DEMAND BALANCE EQUATION (hourly)
******************************************

eq_elc_balance(h,z) ..

* Generation
sum(gen_lim(z,g),var_gen(h,z,g))

* NonVRE Curtailment due to ramp rates
*-sum(non_vre,var_non_vre_curtail(z,h,non_vre))

* Transmission, import-export
-sum(trans_links(z,z_alias,trans),var_trans_flow(h,z_alias,z,trans))
+sum(trans_links(z,z_alias,trans),var_trans_flow(h,z,z_alias,trans)*
    (1-(trans_links_dist_bidir(z,z_alias,trans)*trans_loss(trans))))

$ifThen "%storage%" == ON

* Storage, generated-stored
-sum(s_lim(z,s),var_store(h,z,s))
+sum(s_lim(z,s),var_store_gen(h,z,s))

$endIf

$IF "%pen_gen%" == ON +var_pgen(h,z)

=E= demand(z,h);

* +sum(gen_lim(z,g)$(dem(g)),var_gen(h,z,g));


******************************************

*** Capacity balance ***

eq_new_pcap (g)..
    sum(gen_lim(z,g),var_new_pcap_z(z,g)) =E= var_new_pcap(g);

eq_exist_pcap(g)..
    sum(gen_lim(z,g),var_exist_pcap_z(z,g)) =E= var_exist_pcap(g);

eq_tot_pcap_z(z,g)..
    var_new_pcap_z(z,g) + var_exist_pcap_z(z,g) =E= var_tot_pcap_z(z,g);

eq_tot_pcap(g)..
    sum(z,var_tot_pcap_z(z,g)) =E= var_tot_pcap(g);

*********************
*** VRE equations ***
*********************

* VRE generation is input data x capacity in each region

eq_gen_vre_r(h,vre_lim(vre,z,r))..
    var_vre_gen_r(h,z,vre,r) =L= vre_gen(h,vre,r)*(var_new_vre_pcap_r(z,vre,r)
        +var_exist_vre_pcap_r(z,vre,r));

* VRE gen at regional level aggregated to zonal level

eq_gen_vre(h,z,vre)..
    var_gen(h,z,vre) =E= sum(vre_lim(vre,z,r),var_vre_gen_r(h,z,vre,r));

* VRE capacity across all regions in a zone must be equal to capacity in that
*   zone

eq_new_vre_pcap_z(z,vre)..
    sum(vre_lim(vre,z,r),var_new_vre_pcap_r(z,vre,r)) =E= var_new_pcap_z(z,vre);

eq_exist_vre_pcap_z(z,vre)..
    sum(vre_lim(vre,z,r),var_exist_vre_pcap_r(z,vre,r))
    =E= var_exist_pcap_z(z,vre);

* VRE capacity in each region must be less than or equal to buildable MW as
*   governed by buildable area for each technology in that region

eq_area_max(vre_lim(vre,z,r)) .. var_new_vre_pcap_r(z,vre,r) =L= area(vre,z,r);


*************************
*** NON VRE equations ***
*************************

* Maximum generation of Non VRE

eq_gen_max(gen_lim(z,non_vre),h)$(gen_lin(z,non_vre))..
    var_tot_pcap_z(z,non_vre)*gen_af(non_vre) =G= var_gen(h,z,non_vre);

* Minimum generation of Non VRE

eq_gen_min(mingen_on(z,non_vre),h)$(gen_lin(z,non_vre))..
    var_gen(h,z,non_vre) =G= var_tot_pcap_z(z,non_vre)*gen_mingen_lin(non_vre);

* Ramp equations applied to Non VRE generation, characterised as fraction of
*   total installed capacity per hour

eq_ramp_up(h,ramp_on(z,non_vre))$(gen_lin(z,non_vre))..
    var_gen(h,z,non_vre) =L= var_gen(h-1,z,non_vre)
    +(gen_maxramp(non_vre)*60./gen_unitsize(non_vre))*var_tot_pcap_z(z,non_vre);

eq_ramp_down(h,ramp_on(z,non_vre))$(gen_lin(z,non_vre))..
    var_gen(h,z,non_vre) =G= var_gen(h-1,z,non_vre)
    -(gen_maxramp(non_vre)*60./gen_unitsize(non_vre))*var_tot_pcap_z(z,non_vre);

* Non VRE curtailment due to ramping/min generation

*eq_curtail_max_non_vre(ramp_and_mingen(z,non_vre),h)..
*   var_non_vre_curtail(z,h,non_vre) =L= var_non_vre_gen(z,h,non_vre);


******************************
*** Transmission equations ***
******************************

* Transmitted electricity each hour must not exceed transmission capacity

eq_trans_flow(h,trans_links(z,z_alias,trans))..
    var_trans_flow(h,z,z_alias,trans) =L= var_trans_pcap(z,z_alias,trans);

* Bidirectionality equation is needed when investments into new links are made
*   ...I think :)

eq_trans_bidirect(trans_links(z,z_alias,trans))..
    var_trans_pcap(z,z_alias,trans) =E= var_trans_pcap(z_alias,z,trans);


***********************
*** Misc. equations ***
***********************

* Emissions limit - European average

eq_co2_target(yr)..
    sum((gen_lim(z,non_vre),h)$(hr2yr_map(yr,h)),var_gen(h,z,non_vre)
        *gen_emisfac(non_vre))*1E3
    =L= sum((z,h)$(hr2yr_map(yr,h)),demand(z,h))*%co2intensity%


* Capacity Margin

*scalar dem_max;
*dem_max=smax(h,sum(z,demand(z,h)));

*eq_cap_margin..
*   sum(non_vre,var_tot_pcap(non_vre)*gen_peakaf(non_vre))
*   +sum(vre,var_tot_pcap(vre)*gen_peakaf(vre)) =G= dem_max*1.1;

*eq_max_cap(z,g)..
*   var_cap_z(z,g)+sum(vre_lim(vre,z,r),exist_vre_cap_r(z,vre,r))
*   +gen_exist_pcap_z(z,non_vre) =L= max_cap(z,g)

* Equation for minimum renewable share of generation, set based on restricting
*   non VRE generation.

set flexgen(non_vre) / NaturalgasOCGTnew /;

$IF "%RPS%" == "optimal" $GOTO optimal
    scalar dem_tot;
    dem_tot=sum((z,h),demand(z,h));

    Equations eq_max_non_vre;
    eq_max_non_vre..
        sum((h,z,non_vre)$(gen_lim(z,non_vre) and not flexgen(non_vre)),
            var_gen(h,z,non_vre))
        =E= dem_tot*(1-RPS_scalar);
$label optimal



Model Dispatch /all/;

* don't usually use crossover but can be used to ensure
* a simplex optimal solution is found


Dispatch.OptFile = 1;


$ifThen "%UC%" == ON

Solve Dispatch minimizing costs using MIP;

$else

Solve Dispatch minimizing costs using LP;

$endIf


parameter trans_f(h,z,z_alias,trans);
trans_f(h,z,z_alias,trans)=var_trans_flow.l(h,z_alias,z,trans)
    $(var_trans_flow.l(h,z,z_alias,trans)>1.0);
*display trans_f;

parameter max_bidir_trans;
max_bidir_trans=smax((h,z,z_alias,trans),trans_f(h,z,z_alias,trans));
*display maxtrans;

*parameter pgen_tot;
*pgen_tot=sum((z,h),var_non_vre_gen.l(z,h,"pgen"));

$IF "%log%" == "" $GOTO nolog
scalar now,year,month,day,hour,minute;
now=jnow;
year=gyear(now);
month=gmonth(now);
day=gday(now);
hour=ghour(now);
minute=gminute(now);

scalars  cplex_absgap    absolute gap
         cplex_relgap    relative gap;

cplex_absgap = abs(Dispatch.objest - Dispatch.objval);
cplex_relgap=100*cplex_absgap/abs(Dispatch.objval);

file fname /"%log%"/;
fname.ap=1;
putclose fname "%outname%"","
day:0:0"/"month:0:0"/"year:0:0" "hour:0:0":"minute:0:0","
Dispatch.modelStat:1:0","
Dispatch.resUsd:0
$IF "%UC%" == ON ","cplex_relgap:0
;
$LABEL nolog

* write result parameters

$INCLUDE %codefolderpath%/highres_results.gms

* dump data to GDX
$setEnv GDXCOMPRESS 1
execute_unload "%outname%"

* convert GDX to SQLite


$ifThen "%gdx2sql%" == ON
execute "gdx2sqlite -i %outname%.gdx -o %outname%.db -fast"
$endIf

