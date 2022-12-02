model GlucagonModule
  SinglePool glucagon "in [pg]";
  VariableSource gln_clearance;
  VariableSource gln_release;
  
  parameter Real volume = 11310 "ml";
  parameter Real basal_gln = 72*volume "[pg/ml] * [ml], Sorensen says 80-150";
  parameter Real clearance_k = 910/volume "[ml/min] / [ml]";
  parameter Real basal_release = clearance_k * basal_gln "TODO pg/min";
  
  // effects of glucose and insulin, inputs
  Real insulin;
  Real glucose;
  
  SlowSaturation ins_delay(k=200*60);

  
  Real ins_eff = 1.31 - 0.61*tanh(1.06*(ins_delay.out-0.47));
  Real glu_eff = 2.93 - 2.10*tanh(4.18*(glucose-0.61));
    
  // output variables
  Real gln_normed = glucagon.metabolite_mass/basal_gln;
  Real glucagon_C = glucagon.metabolite_mass/volume;
  
initial equation
  glucagon.metabolite_mass = basal_gln;
  
equation
  ins_delay.inp = insulin;

  // clearance
  gln_clearance.k_source = -clearance_k * glucagon.metabolite_mass;
  
  connect(glucagon.metabolite, gln_clearance.metabolite);


  // release
  gln_release.k_source = glu_eff * ins_eff * basal_release;
  
  connect(glucagon.metabolite, gln_release.metabolite);

end GlucagonModule;
