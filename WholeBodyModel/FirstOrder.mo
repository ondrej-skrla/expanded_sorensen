model FirstOrder
  extends SinkToSource;
    
  parameter Real k "reaction rate [mg/min]";
  parameter Real basal_inp "basal mass of input [mg]";

protected
  Real inp_normed = met_inp.M/basal_inp "normalized input";
  
equation 
  60*met_inp.Q = k*inp_normed;

end FirstOrder;
