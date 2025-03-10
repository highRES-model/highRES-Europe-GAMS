******************************************
* highres unit commitment module
******************************************

option IntVarUp=0
$ONEPS

* f_response = 10 seconds
* reserve = 20 minutes



* vre forecast errors (fraction of output)
* - both from NatGrid 2017 - Demand Forecasting
parameter
vre_margin(vre)
/
Windonshore 0.14
Windoffshore 0.14
Solar 0.12
/
;





set map_minup(h,non_vre,h)
    map_mindown(h,non_vre,h);


*map1_minup(h,non_vre,h_alias) = ord(h_alias) ge (ord(h) - gen_minup(non_vre)+1)
*and ord(h_alias) lt ord(h);
map_minup(h,non_vre,h_alias+(ord(h)-gen_minup(non_vre)))$[hh_minup(h_alias) and
    ord (h_alias)<gen_minup(non_vre)] = yes;

*map1_mindown(h,non_vre,h_alias) = ord(h_alias) ge (ord(h)
*    - gen_mindown(non_vre)+1) and ord(h_alias) lt ord(h);
map_mindown(h,non_vre,h_alias+(ord(h)-gen_mindown(non_vre)))
    $[hh_mindown(h_alias) and ord (h_alias)<gen_mindown(non_vre)] = yes;

$ontext

set diff1(h,non_vre,h)
diff2(h,non_vre,h);

diff1(h,non_vre,h_alias) = map1_minup(h,non_vre,h_alias) xor
    map_minup(h,non_vre,h_alias);
diff2(h,non_vre,h_alias) = map1_mindown(h,non_vre,h_alias) xor  
    map_mindown(h,non_vre,h_alias);
* abort$(card(diff))

* 'map1 and map2 not identical'
display diff1,diff2;

$offtext

parameter gen_max_res(non_vre,service_type);

* compute maximum power ramp in MW for each tech for in each reserve window

gen_max_res(non_vre,"reserve")$(gen_uc_int(non_vre) or gen_uc_lin(non_vre)) = 
    gen_maxramp(non_vre)*res_time;
gen_max_res(non_vre,"f_response")$(gen_uc_int(non_vre) or gen_uc_lin(non_vre)) =
    gen_maxramp(non_vre)*f_res_time;

$ontext
max_res("Nuclear","f_response")=1.7;
max_res("NaturalgasCCGTwithCCSnewOT","f_response")=8.;
max_res("NaturalgasOCGTnew","f_response")=14.5;

* CCGT = 47.0 MW/min (=> (47/2)/500 = 0.05 rounded up), OCGT = 87 MW/min
* taken from https://www.sciencedirect.com/science/article/pii/S1364032117309206
$offtext

* rescale from MW to GW for better numerics (allegedly)
* also need to change var_H equation to include a /1E3 and remove
* var_freq_req scaling

* startupcost was in £k, now in £m as objective function is now in £m


gen_startupcost(non_vre)=gen_startupcost(non_vre)/MWtoGW;
unit_cap_lim_z=unit_cap_lim_z/MWtoGW;

parameter res_req(h)                     operating reserve requirement based on
                                         demand level;
res_req(h)=sum(z$(uc_z(z)),demand(z,h))*res_margin;


*************************

*************************

Positive variables
var_tot_n_units_lin(z,non_vre)              total number of units (linear)
var_new_n_units_lin(z,non_vre)              number of new units (linear)
var_exist_n_units_lin(z,non_vre)            number of existing units (linear)
var_up_units_lin(h,z,non_vre)               units starting up by tech zone hour (linear)
var_down_units_lin(h,z,non_vre)             units shutdown by tech zone hour (linear)
var_com_units_lin(h,z,non_vre)              units committed by tech zone hour (linear)
var_res(h,z,non_vre)                        operating reserve offered by tech zone hour (MW)
var_res_quick(h,z,non_vre)                  quick start operating reserve by tech zone hour (MW)
$IF "%f_res%" == ON var_f_res(h,z,non_vre)  frequency resonse by tech zone hour (MW)
;

Integer variables
var_tot_n_units(z,non_vre)                  total units per zone
var_new_n_units(z,non_vre)                  newly installed units per zone
var_exist_n_units(z,non_vre)                existing units per zone
var_up_units(h,z,non_vre)                   units started per tech hour and zone
var_down_units(h,z,non_vre)                 units shutdown per tech hour and zone
var_com_units(h,z,non_vre)                  committed units per tech hour and zone
;





* set integer upper limit for existing capacity being represented as units,
* based on:
* a) existing cap, floored to ensure it doesn't breach existing cap limit
* b) if no existing cap, set equal to 0

var_exist_n_units.UP(z,non_vre)$(gen_uc_int(non_vre) and
    gen_exist_pcap_z(z,non_vre,"UP") and uc_z(z))
    =floor(gen_exist_pcap_z(z,non_vre,"UP")/gen_unitsize(non_vre));
var_exist_n_units.L(z,non_vre)$(gen_uc_int(non_vre) and 
    gen_exist_pcap_z(z,non_vre,"UP") and uc_z(z))
    =floor(gen_exist_pcap_z(z,non_vre,"UP")/gen_unitsize(non_vre));

var_exist_n_units.LO(z,non_vre)$(gen_uc_int(non_vre) and
    gen_exist_pcap_z(z,non_vre,"LO") and uc_z(z))
    = floor(gen_exist_pcap_z(z,non_vre,"LO")/gen_unitsize(non_vre));
var_exist_n_units.FX(z,non_vre)$(gen_uc_int(non_vre) and
    gen_exist_pcap_z(z,non_vre,"FX") and uc_z(z))
    = floor(gen_exist_pcap_z(z,non_vre,"FX")/gen_unitsize(non_vre));

var_exist_n_units.UP(z,non_vre)$(gen_uc_int(non_vre) and not
    sum(lt,gen_exist_pcap_z(z,non_vre,lt)) and uc_z(z))=0.;

* set integer upper limit for new capacity being represented as units, based on:
* a) total cap limit per zone
* b) some arbitrary number (unit_cap_lim_z)

var_tot_n_units.UP(z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    gen_lim_pcap_z(z,non_vre,'UP') < INF and uc_z(z) )
    =ceil(gen_lim_pcap_z(z,non_vre,'UP')/gen_unitsize(non_vre));
var_tot_n_units.UP(z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    gen_lim_pcap_z(z,non_vre,'UP')
    = INF and uc_z(z))=ceil(unit_cap_lim_z/gen_unitsize(non_vre));

*var_n_units.UP(z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre))
*=ceil(50000/gen_unitsize(non_vre));
var_up_units.UP(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))=var_tot_n_units.UP(z,non_vre);
var_down_units.UP(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))=var_tot_n_units.UP(z,non_vre);
var_com_units.UP(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))=var_tot_n_units.UP(z,non_vre);

var_com_units.UP(h,z,non_vre)$(gen_uc_int(non_vre) and not gen_lim(z,non_vre)
    and uc_z(z))=0.;

Equations

* integer operability equations

eq_uc_tot_units
eq_uc_units
eq_uc_unit_state
eq_uc_cap
eq_uc_exist_cap
eq_uc_gen_max
eq_uc_gen_min
eq_uc_gen_minup
eq_uc_gen_mindown

* linear operability equations

eq_uc_tot_units_lin
eq_uc_units_lin
eq_uc_unit_state_lin
eq_uc_cap_lin
eq_uc_exist_cap_lin
eq_uc_gen_max_lin
eq_uc_gen_min_lin
eq_uc_gen_minup_lin
eq_uc_gen_mindown_lin

* reserves/response

$ifThen "%f_res%" == ON

eq_uc_response
eq_uc_max_response
eq_uc_max_response_lin

$endIf

eq_uc_reserve_quickstart
eq_uc_reserve


* no quickstart frequency response at the moment
*eq_uc_response_quickstart

eq_uc_max_reserve
eq_uc_max_reserve_lin



;

** Capacity balance equations

eq_uc_tot_units(z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_tot_n_units(z,non_vre)
        =E= var_exist_n_units(z,non_vre)+var_new_n_units(z,non_vre);

eq_uc_cap(z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and uc_z(z))..
    var_new_pcap_z(z,non_vre)
    =E= var_new_n_units(z,non_vre)*gen_unitsize(non_vre);

eq_uc_exist_cap(z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_exist_pcap_z(z,non_vre)
        =E= var_exist_n_units(z,non_vre)*gen_unitsize(non_vre);


eq_uc_tot_units_lin(z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_tot_n_units_lin(z,non_vre)
        =E= var_exist_n_units_lin(z,non_vre)+var_new_n_units_lin(z,non_vre);

eq_uc_cap_lin(z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_new_pcap_z(z,non_vre)
        =E= var_new_n_units_lin(z,non_vre)*gen_unitsize(non_vre);

eq_uc_exist_cap_lin(z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_exist_pcap_z(z,non_vre)
        =E= var_exist_n_units_lin(z,non_vre)*gen_unitsize(non_vre);

** Integer operability

* total committed units must be less than installed units

eq_uc_units(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_com_units(h,z,non_vre)
        =L= var_tot_n_units(z,non_vre);

* committment state of units

eq_uc_unit_state(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_com_units(h,z,non_vre) =E= var_com_units(h-1,z,non_vre)
            +var_up_units(h,z,non_vre)-var_down_units(h,z,non_vre);

* maximum generation limit

eq_uc_gen_max(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_com_units(h,z,non_vre)*gen_unitsize(non_vre)*gen_af(non_vre)
        =G= var_gen(h,z,non_vre) + var_res(h,z,non_vre)
$IF "%f_res%" == ON +var_f_res(h,z,non_vre)
        ;

* minimum stable generation

eq_uc_gen_min(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_gen(h,z,non_vre) =G= var_com_units(h,z,non_vre)*gen_mingen(non_vre)
            *gen_unitsize(non_vre);

* minimum up time

eq_uc_gen_minup(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    (gen_minup(non_vre) > 1) and (ord(h) > 1) and uc_z(z))..
        sum(map_minup(h,non_vre,h_alias),var_up_units(h_alias,z,non_vre))
        =L= var_com_units(h,z,non_vre);

* minimum down time

eq_uc_gen_mindown(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    (gen_mindown(non_vre) > 1) and (ord(h) > 1) and uc_z(z))..
        sum(map_mindown(h,non_vre,h_alias),var_down_units(h_alias,z,non_vre))
        =L= var_tot_n_units(z,non_vre)-var_com_units(h,z,non_vre);


** Linear operability


eq_uc_units_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_com_units_lin(h,z,non_vre) =L= var_tot_n_units_lin(z,non_vre);

eq_uc_unit_state_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)
    and uc_z(z))..
        var_com_units_lin(h,z,non_vre) =E= var_com_units_lin(h-1,z,non_vre)
            +var_up_units_lin(h,z,non_vre)-var_down_units_lin(h,z,non_vre);

eq_uc_gen_max_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)
    and uc_z(z))..
        var_com_units_lin(h,z,non_vre)*gen_unitsize(non_vre)*gen_af(non_vre)
        =G= var_gen(h,z,non_vre) + var_res(h,z,non_vre)
$IF "%f_res%" == ON +var_f_res(h,z,non_vre)
        ;

eq_uc_gen_min_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_gen(h,z,non_vre) =G= var_com_units_lin(h,z,non_vre)*
            gen_mingen(non_vre)*gen_unitsize(non_vre)*gen_af(non_vre);

eq_uc_gen_minup_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre) and
    (gen_minup(non_vre) > 1) and (ord(h) > 1) and uc_z(z))..
        sum(map_minup(h,non_vre,h_alias),var_up_units_lin(h_alias,z,non_vre))
        =L= var_com_units_lin(h,z,non_vre);

eq_uc_gen_mindown_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)
    and (gen_mindown(non_vre) > 1) and (ord(h) > 1) and uc_z(z))..
        sum(map_mindown(h,non_vre,h_alias),
            var_down_units_lin(h_alias,z,non_vre))
        =L= var_tot_n_units_lin(z,non_vre)-var_com_units_lin(h,z,non_vre);



** Reserves

* Quickstart - units which can come online and ramp to full in the reserve
*   window -> OCGT only

eq_uc_reserve_quickstart(h,z,non_vre)$(gen_lim(z,non_vre) and gen_quick(non_vre)
    and uc_z(z))..
        (var_tot_n_units_lin(z,non_vre)-var_com_units_lin(h,z,non_vre))*
            gen_unitsize(non_vre)*gen_af(non_vre)
        =G= var_res_quick(h,z,non_vre);

* Max reserve potential if needed - can be used to simulate different time
*   scales over which reserve is offered

eq_uc_max_reserve(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_res(h,z,non_vre) =L= var_com_units(h,z,non_vre)
            *gen_max_res(non_vre,"reserve")*gen_af(non_vre);

eq_uc_max_reserve_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)
    and not gen_quick(non_vre) and uc_z(z))..
        var_res(h,z,non_vre) =L= var_com_units_lin(h,z,non_vre)
            *gen_unitsize(non_vre)*gen_af(non_vre)
            *gen_max_res(non_vre,"reserve");

$ifThen "%f_res%" == ON

* Max frequency response potential

eq_uc_max_response(h,z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre) and
    uc_z(z))..
        var_f_res(h,z,non_vre)=L= var_com_units(h,z,non_vre)
            *gen_max_res(non_vre,"f_response")*gen_af(non_vre);

eq_uc_max_response_lin(h,z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)
    and uc_z(z))..
        var_f_res(h,z,non_vre) =L= var_com_units_lin(h,z,non_vre)
            *gen_max_res(non_vre,"f_response")*gen_af(non_vre);

$endIf


$ifThen.a "%f_res%" == ON

equations
eq_uc_H
;

Positive variables
var_H(h)                                 system inertia per hour (GWs per Hz)
;


scalar f_0 / 50. /;
scalar p_loss / 1650. /;
scalar p_loss_inertia  / 7./;

*positive variable var_inertia(h);

*+var_inertia(h)

* compute system inertia

eq_uc_H(h)..
    var_H(h) =E= sum((z,non_vre)$(gen_uc_int(non_vre) and gen_lim(z,non_vre)
                    and uc_z(z)),gen_inertia(non_vre)*var_com_units(h,z,non_vre)
                    *(gen_unitsize(non_vre)*MWtoGW/1E3)/f_0)
                    +sum((z,non_vre)$(gen_uc_lin(non_vre) and gen_lim(z,non_vre)
                        and uc_z(z)),gen_inertia(non_vre)
                        *var_com_units_lin(h,z,non_vre)
                        *(gen_unitsize(non_vre)*MWtoGW/1E3)/f_0)
$ifthen.b "%store_uc%" == ON
    +sum((z,s)$(s_lim(z,s) and uc_z(z) and store_uc_lin(s)),
    store_inertia(s)*var_store_com_units_lin(h,z,s)*(store_unitsize(s)*MWtoGW/1E3)/f_0);
$endif.b    

* -((p_loss/1E3)*p_loss_inertia/f_0);

var_H.LO(h) = 0.825;
*var_H.LO(h) = 0.825*2;
*var_H.UP(h) = 9..


set seg /1*8/;
set lin_param /slope,intercept/;

table linearise(seg,lin_param)
        slope        intercept
1        -3.4        10.21
2        -1.13        5.67
3        -0.57        3.97
4        -0.34        3.06
5        -0.23        2.5
6        -0.16        2.11
7        -0.12        1.82
8        -0.09        1.61


;


*var_H.UP(h) = 9.;
*var_H.LO(h) = 1.;



equations
eq_uc_freq_req
;

positive variable var_freq_req(h);

* Frequency nadir has to be reached by delivery time => 5 seconds in this case
* issue with setting delivery time to 1 second means Td < Tn, below ensures that
* Td >= Tn

var_freq_req.LO(h)=p_loss/MWtoGW;

*var_freq_req.UP(h)=var_H.LO(h)*linearise("1","slope")
*   +linearise("1","intercept");
*var_freq_req.LO(h)=var_H.UP(h)*linearise("6","slope")
*   +linearise("6","intercept");

* piecewise linear approximation of frequency response requirement.
* see:
* Teng, F. et al. (2017). Full stochastic scheduling for low-carbon electricity
*   systems.
*   IEEE Transactions on Automation Science and Engineering, 14(2), 461-470

eq_uc_freq_req(h,seg)..
    var_freq_req(h) =G= (linearise(seg,"slope")*var_H(h)
        +linearise(seg,"intercept"))*1E3/MWtoGW;


$endIf.a

$ifThen.a "%storage%" == ON


Equations
eq_uc_store_res_level
eq_uc_store_res_max
$IF "%f_res%" == ON eq_uc_store_f_res_max
;

* limits response and reserve offered by storage to be less than the storage
*   level - only for techs which can offer response or reserve

eq_uc_store_res_level(h,s_lim(z,s))$((store_max_res(s) > 0. or
    store_max_freq(s) > 0.) and uc_z(z))..
        var_store_res(h,z,s)
$IF "%f_res%" == ON +var_store_f_res(h,z,s)
        =L= var_store_level(h,z,s);

* limits reserve offered by storage based on how much of its capacity can come
*   online within the reserve time window - only for techs which can offer
*   reserve

* store_max_res(s) > 0. and

eq_uc_store_res_max(h,s_lim(z,s))$(uc_z(z) and not store_uc_lin(s))..
    var_store_res(h,z,s)
    =L= var_tot_store_pcap_z(z,s)*store_af(s)*store_max_res(s);

* limits resposne offered by storage based on how much of its capacity can come
*   online within the response window - only for techs which can offer response

* store_max_freq(s) > 0. and

$ifThen.b "%f_res%" == ON

eq_uc_store_f_res_max(h,s_lim(z,s))$(uc_z(z) and not
store_uc_lin(s))..
var_store_f_res(h,z,s) 
=L= var_tot_store_pcap_z(z,s)*store_af(s)*store_max_freq(s);

$endif.b

$endIf.a


* Main reserve equation

eq_uc_reserve(h) ..
*   spinning - only techs which can offer reserve
    sum((z,non_vre)$(gen_lim(z,non_vre) and gen_max_res(non_vre,"reserve") > 0.
    and uc_z(z)),var_res(h,z,non_vre))

*   quick start
    +sum((z,non_vre)$(gen_lim(z,non_vre) and gen_quick(non_vre) and uc_z(z)),
        var_res_quick(h,z,non_vre))
        
$ifthen "%store_uc%" == ON       
    +sum((z,s)$(s_lim(z,s) and store_uc_lin(s) and
      uc_z(z) and store_quick(s)),
      var_store_res_quick(h,z,s))
$endif


*   storage

$ifthen "%storage%" == ON
    +sum(s_lim(z,s)$(store_max_res(s) > 0. and uc_z(z) 
      and not store_uc_lin(s)),var_store_res(h,z,s))
    +sum(s_lim(z,s)$(uc_z(z) and store_uc_lin(s)),
      var_store_res(h,z,s))
$endif      

    -sum((z,vre)$(uc_z(z)),var_gen(h,z,vre)*vre_margin(vre))

    =G= res_req(h);




$ifThen.a "%f_res%" == ON

* Main response equation
eq_uc_response(h)..

*   spinning
        sum((z,non_vre)$(gen_lim(z,non_vre) and
            gen_max_res(non_vre,"f_response") > 0. and uc_z(z)),
            var_f_res(h,z,non_vre))

*   storage
$ifthen.b "%storage%" == ON
    +sum(s_lim(z,s)$(store_max_freq(s) > 0. and
    uc_z(z) and not store_uc_lin(s)),
    var_store_f_res(h,z,s))
    +sum(s_lim(z,s)$(uc_z(z) and store_uc_lin(s)),
    var_store_f_res(h,z,s))
$endif.b

        =G= var_freq_req(h);

$endIf.a
