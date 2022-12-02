model InputMassFlow "Input metabolite mass flow from inflowing organs"
  extends VectorSource;
  
  RealInput vessel_in_C[n_inp, MPlasma] "concentrations from input organs [mg/dl]";
  parameter Real in_blood_flows[n_inp, MPlasma] "blood flows from input organs [dl/min]";
  
  Real in_flow_rate[MPlasma] = fill(1, n_inp) * (in_blood_flows .* vessel_in_C) "input flow rate [mg/min]";

equation
  k_source = in_flow_rate;

end InputMassFlow;
