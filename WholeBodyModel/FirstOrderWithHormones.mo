model FirstOrderWithRegulation
  extends SinkToSource;
    
  parameter Real k "reaction rate [mg/min]";
  parameter Real basal_inp "basal mass of input [mg]";
  Real reg_eff "the effect of regulation [-]";

protected
  Real inp_normed = met_inp.M/basal_inp "normalized input";
  
equation 
  60*met_inp.Q = reg_eff*k*inp_normed;
  
end FirstOrderWithRegulation;
