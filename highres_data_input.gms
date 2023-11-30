$ONEPS
$ONEMPTY

Sets

lt / UP, LO, FX /


r regions /
$BATINCLUDE %datafolderpath%/%vre_restrict%_regions.dd
/

z zones /
$BATINCLUDE %datafolderpath%/zones.dd
/
;


$INCLUDE %datafolderpath%/%weather_yr%_temporal.dd

alias(h,h_alias);

alias(z,z_alias) ;


* these sets have to be read in here, they can't be read in at the same time as
* parameter definitions

set g;
set non_vre(g);
set vre(g);
set trans;

parameter gen_capex(g);
parameter gen_varom(g);
parameter gen_fom(g);
parameter gen_fuelcost(g);
parameter gen_mingen(g);
parameter gen_emisfac(g);
parameter gen_maxramp(g);
parameter gen_af(g);
parameter gen_peakaf(g);
parameter gen_cap2area(g);
parameter gen_lim_pcap_z(z,g,lt);
parameter gen_lim_ecap_z(z,g,lt);
parameter gen_exist_pcap_z(z,g,lt);
parameter gen_exist_ecap_z(z,g,lt);
parameter gen_fx_natcap(g);

parameter gen_unitsize(non_vre);
parameter gen_startupcost(non_vre);
parameter gen_minup(non_vre);
parameter gen_mindown(non_vre);
parameter gen_inertia(non_vre);

parameter trans_links_cap(z,z_alias,trans);
parameter trans_links_dist(z,z_alias,trans);
parameter trans_loss(trans);
parameter trans_varom(trans);
parameter trans_line_capex(trans);
parameter trans_sub_capex(trans);

parameter area(vre,z,r);
parameter vre_gen(h,vre,r);

parameter demand(z,h);

scalar co2_budget;


$INCLUDE %datafolderpath%/%psys_scen%_gen.dd
$INCLUDE %datafolderpath%/trans.dd
$INCLUDE %datafolderpath%/%esys_scen%_co2_budget.dd

* need to switch between agg and not for areas currently

$INCLUDE %datafolderpath%/vre_areas_%weather_yr%_%vre_restrict%.dd
*$BATINCLUDE %datafolderpath%/vre_%weather_yr%_agg_new_%area_scen%.dd
$INCLUDE %datafolderpath%/%esys_scen%_demand_%dem_yr%.dd

$gdxin %datafolderpath%/vre_%weather_yr%_%vre_restrict%.gdx
$load vre_gen

* $IF %fx_natcap% == "YES" $INCLUDE %datafolderpath%/%esys_scen%_gen_fx_natcap.dd

*gen_capex(g)=gen_capex%model_yr%(g);
*gen_fuelcost(g)=gen_fuelcost%model_yr%(g);

parameter exist_vre_cap_r(vre,z,r);

*$BATINCLUDE %datafolderpath%/vre_per_zone_2016.dd

