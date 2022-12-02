model GLYtoG6P "glycogen to glucose-6-phosphate reaction"
  extends SinkToSource;
    
  parameter Real kc_7 "Maximum speed of reaction [mg/min]";
  parameter Real kc_8 "Michealis-Menten constant [-]";
  parameter Real basal_inp "basal mass of input [mg]";
  Real ins_eff "the effect of insulin [-]";

protected
  Real inp_normed = met_inp.M/basal_inp "normalized input";
  
equation        
  60*met_inp.Q = ins_eff*(kc_7*inp_normed/(kc_8 + inp_normed));
  
end GLYtoG6P;
