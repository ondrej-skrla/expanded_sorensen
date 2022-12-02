model OutputMassFlow "organ output metabolite mass flow"
  extends VectorSource;
  
  parameter Real out_blood_flow[MPlasma] "output blood flows for individual metabolites [dl/min]";
  
  RealInput vessel_out_C[MPlasma] "output concentration of metabolites [mg/dl]";
  Real out_flow_rate[MPlasma] = out_blood_flow .* vessel_out_C "[mg/min]";

equation
  k_source = -out_flow_rate;
  
end OutputMassFlow;
