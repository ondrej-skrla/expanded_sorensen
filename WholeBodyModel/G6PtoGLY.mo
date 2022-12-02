model G6PtoGLY "glucose-6-phosphate to glycogen reaction"
  extends SinkToSource;
   
  parameter Real kc_4 "basal glycogenesis rate [mg/min]";
  parameter Real kc_5 "maximum glycogen storage [mg]";
  parameter Real kc_6 "tanh steepness of slope [mg]";
  parameter Real basal_inp "basal mass of input [mg]";
  Real reg_eff "the effect of regulation [-]";
  
protected
  Real inp_normed = met_inp.M/basal_inp "normalized input [-]";

equation
        
  60*met_inp.Q = reg_eff*kc_4*inp_normed*0.5*(1+tanh(kc_6*(kc_5-met_out.M)));

end G6PtoGLY;
