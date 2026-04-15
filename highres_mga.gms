******************************************
* highres MGA module
******************************************

$ONEPS
$ONEMPTY

* Parameters such as slack and cost-optimal value.
$INCLUDE %mgapath%/mga_parameters.dd

* Subsets for focus of MGA objective
Set
    mga_z(z)  "zones included in MGA objective"
    mga_g(g)  "technologies included in MGA objective"
;

mga_z(z) = yes$(par_mga_z(z));
mga_g(g) = yes$(par_mga_g(g));

Variables
mga_objective
;

Equations
eq_mga_obj
eq_max_costs
;

* Sum over the subsets 
eq_mga_obj..
    mga_objective =E=
        sum( (mga_z(z), mga_g(g)),
              var_exist_pcap_z(z,g)
            + var_new_pcap_z(z,g)
        );

eq_max_costs..
    costs =L= par_optimal_cost * (1+par_slack);
