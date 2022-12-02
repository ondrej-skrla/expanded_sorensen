model AdiposeTissue

  // specify the base model
  extends OrganModel
    (n_inp=1,
     separate=true,
     extra_cell_volume = AdiposeConstants.volume_non_cellular,
     intra_cell_volume = AdiposeConstants.volume_cellular,
                     
     input_mass_flow.in_blood_flows = {AdiposeConstants.input_blood_flow},
     output_mass_flow.out_blood_flow = AdiposeConstants.input_blood_flow,
                     

     membrane_transport.extra_volume = AdiposeConstants.volume_non_cell_exchange,

     extra_cell_C(start=AdiposeConstants.basal_C_out, each fixed=true),
     transfer_C(start=AdiposeConstants.basal_C_transfer, each fixed=true),       
     locked_C(start=AdiposeConstants.basal_C_locked, each fixed=true)
     );
     
  // basal masses
  parameter Real[MChange] basal_trans = AdiposeConstants.basal_C_transfer * intra_cell_volume "[mg]";
  parameter Real[MLocked] basal_locked = AdiposeConstants.basal_C_locked * intra_cell_volume "[mg]";
            
  // Glucose reactions
  GLUtoG6P glu_to_g6p(kc_1 = 100,
                      kc_2 = 9,
                      kc_3 = 20,
                      basal_inp = basal_trans[MChange.GLU],
                      basal_out = basal_locked[MLocked.G6P]);
                                          
  FirstOrderWithRegulation
             g6p_to_pyr(k = 10,
                        basal_inp = basal_locked[MLocked.G6P]);
  
  FirstOrder pyr_to_lac(k = 10,
                        basal_inp = basal_locked[MLocked.PYR]);
                        
  FirstOrder lac_to_pyr(k = 5,
                        basal_inp = basal_trans[MChange.LAC]);
  
  TanhInhibWithRegulation
             pyr_to_acoa(k_1 = 5,
                         k_2 = 1,
                         k_3 = 1,
                         basal_inp = basal_locked[MLocked.PYR],
                         basal_out = basal_locked[MLocked.ACOA]);
                               
  
  ACOA_SINK
             acoa_sink(k_1 = 5,
                       k_2 = 1.5,
                       k_3 = 0.8,
                       basal_inp = basal_locked[MLocked.ACOA]);
                         
  // Glycerol source from lipolysis
  VariableSource glycerol_cell_source;
  parameter Real glc_basal_release = 9 "[mg/min]"; 
  
  // effects
  
  VariableSource ins_transfer "transfer from blood to interstitium";
  
  parameter Real volume_I = 6.74 "dl";
  parameter Real k_trans_ins = 20 "min";
  parameter Real basal_ins_C = 0.279 "mU/dl";
  Real insulin_venous = extra_cell_C[MPlasma.INS] "venous insulin";
  Real insulin "interstitial insulin";
  Real insulin_normed = insulin/basal_ins_C;
  
  // applies to membrane transport, glu to g6p and g6p to pyr
  Real glu_uptake_ins_eff = 2.3+1.4*tanh(0.338*(insulin-5.82));
  
  // applies to pyr to acoa
  Real oxidation_ins_eff = 1.74+0.8*tanh(0.338*(insulin-5.82));
  
  // effect of insulin on lipolysis and therefore glycerol release
  Real glc_ins_eff = 1.1+0.7*tanh(0.8*(0.8-insulin));
  

  
initial equation

  insulin = basal_ins_C;
  
equation
  
  
  // Insulin membrane transfer
  
  ins_transfer.k_source = -volume_I*(insulin_venous-insulin)/k_trans_ins;
  
  60*volume_I*der(insulin) = volume_I*(insulin_venous-insulin)/k_trans_ins -
                             1.37*insulin;  
                             
  connect(ins_transfer.metabolite, extra_pool.metabolites[MPlasma.INS]);

  // Reactions occuring in cells and interstitium_pool
  
  // Insulin effect on membrane transport (GLUT4 translocation)
  membrane_transport.k_transport = AdiposeConstants.transport_const .*
                                   {glu_uptake_ins_eff,
                                   1, 1, 1};
  
  // Glycerol source from lipolysis
  connect(transfer_pool.metabolites[MChange.GLC], glycerol_cell_source.metabolite);
  glycerol_cell_source.k_source = glc_basal_release * glc_ins_eff;
  
  // Glucose reactions
  
  // GLU to G6P
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p.met_out);
  glu_to_g6p.reg_eff = glu_uptake_ins_eff;
  
  // G6P to PYR
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], g6p_to_pyr.met_out);
  g6p_to_pyr.reg_eff = glu_uptake_ins_eff;
  
  // PYR to LAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], pyr_to_lac.met_out);
  
  // LAC to PYR
  connect(transfer_pool.metabolites[MChange.LAC], lac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], lac_to_pyr.met_out);
  
  // PYR to ACoA
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_acoa.met_inp);
  connect(locked_pool.metabolites[MLocked.ACOA], pyr_to_acoa.met_out);
  pyr_to_acoa.reg_eff = oxidation_ins_eff;
  
  // ACoA to energy (sink)
  connect(locked_pool.metabolites[MLocked.ACOA], acoa_sink.met_inp);

end AdiposeTissue;
