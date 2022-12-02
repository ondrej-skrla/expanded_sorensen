model LiverTissue

  // specify the base model
  extends OrganModel
    (n_inp=2,
     separate=true,
     extra_cell_volume = LiverConstants.volume_non_cellular,
     intra_cell_volume = LiverConstants.volume_cellular,
                     
     input_mass_flow.in_blood_flows = LiverConstants.input_blood_flow,
     output_mass_flow.out_blood_flow = LiverConstants.input_blood_flow[1] +
                                       LiverConstants.input_blood_flow[2],
                     

     membrane_transport.extra_volume = LiverConstants.volume_non_cell_exchange,
     membrane_transport.k_transport = LiverConstants.transport_const,
                     
     //hormone_mass(start=start_hormones, each fixed=true),
     //basal_hormones = LiverConstants.basal_hormones,
     
     extra_cell_C(start=LiverConstants.basal_C_out, each fixed=true),
     transfer_C(start=LiverConstants.basal_C_transfer, each fixed=true),       
     locked_C(start=LiverConstants.basal_C_locked, each fixed=true)
     );
       
  // basal masses
  parameter Real[MChange] basal_trans = LiverConstants.basal_C_transfer * intra_cell_volume "[mg]";
  parameter Real[MLocked] basal_locked = LiverConstants.basal_C_locked * intra_cell_volume "[mg]";
  
  // Insulin release from pancreas
  InsulinRelease ins_release;
  VariableSource pancreas_insulin_inflow;
  
  // Insulin clearance
  VariableSource ins_clearance;
    
                     
  // Glycolysis reactions
  // GLU to G6P, hexokinase
  GLUtoG6P glu_to_g6p(kc_1 = 27,
                      kc_2 = 0.0125,
                      kc_3 = 1,
                      basal_inp = basal_trans[MChange.GLU],
                      basal_out = basal_locked[MLocked.G6P]);
                      
  // GLU to G6P, glucokinase
  MichealisMentenWithReg
             glu_to_g6p_glucokinase
                       (v_max = 325,
                        k_m = 1.5,
                        basal_inp = basal_trans[MChange.GLU]);
                                          
  //G6PtoG3P
  FirstOrder 
             g6p_to_g3p(k = 20,
                        basal_inp = basal_locked[MLocked.G6P]);
  
  //G3PtoPYR
  FirstOrderWithRegulation
             g3p_to_pyr(k = 20,
                        basal_inp = basal_locked[MLocked.G3P]);
                        
  //PYRtoLAC
  FirstOrder pyr_to_lac(k = 10,
                        basal_inp = basal_locked[MLocked.PYR]);
  
  //PYRtoACOA
  TanhInhibWithRegulation
             pyr_to_acoa(k_1 = 20,
                         k_2 = 0.1,
                         k_3 = 1,
                         basal_inp = basal_locked[MLocked.PYR],
                         basal_out = basal_locked[MLocked.ACOA]);
  
  // ACOA_SINK
  ACOA_SINK acoa_burn(k_1 = 20,
                      k_2 = 2.5,
                      k_3 = 0.43,
                      basal_inp = basal_locked[MLocked.ACOA]);
                      
  
  // Glycogen cycling
  // G6P to GLY
  TanhInhibWithRegulation
           g6p_to_gly(k_1 = 27*0.5,
                      k_2 = 2,
                      k_3 = 10,
                      basal_inp = basal_locked[MLocked.G6P],
                      basal_out = basal_locked[MLocked.GLY]);
  
  // GLY to G6P
  TanhWithRegulation
             gly_to_g6p(k_1 = 250,
                        k_2 = 10,
                        basal_inp = basal_locked[MLocked.GLY]);  
  // Gluconeogenesis reactions
  // AACtoPYR
  FirstOrder 
             aac_to_pyr(k = 17,
                        basal_inp = basal_trans[MChange.AAC]);
  
  // LACtoPYR
  FirstOrder
             lac_to_pyr(k = 37,
                        basal_inp = basal_trans[MChange.LAC]);
  
  // PYRtoG3P
  FirstOrder
             pyr_to_g3p(k = 44,
                        basal_inp = basal_locked[MLocked.PYR]);
  
  // GLCtoG3P
  FirstOrder
             glc_to_g3p(k = 6,
                        basal_inp = basal_trans[MChange.GLC]);
  
  // G3PtoG6P
  FirstOrder
             g3p_to_g6p(k = 50,
                        basal_inp = basal_locked[MLocked.G3P]);
  
  // G6PtoGLU
  MichealisMentenWithReg
             g6p_to_glu(v_max = 945,
                        k_m = 2.5,
                        basal_inp = basal_locked[MLocked.G6P]);

  // normalized glucose, insulin for regulation
  Real glu_normed = transfer_pool.metabolites[MChange.GLU].M / basal_trans[MChange.GLU];
  // Real ins_normed = extra_pool.metabolites[MPlasma.INS].M / 
  
  Real locked_normed[MLocked] = locked_pool.metabolites.M ./ basal_locked;

  
  
  // effects
  Real insulin = extra_cell_out_C[MPlasma.INS] / 
    LiverConstants.basal_C_out[MPlasma.INS];
  Real glucagon;
  
  // function for peak decline in glycogenolysis
  PeakDecline gln_glycolysis(k_dec = 65*60,
                            stable_lev = 0.4);
                            
  
  // function for slow rise in glucagon gluconeogenesis
  SlowSaturation gln_gluconeo(k = 40*60);     
  
  // slow rise in insulin effect
  SlowSaturation ins_effect(k = 20*60);                  

  
  // GLU to G6P glucokinase
  Real glu_to_g6p_ins_eff = 0.77*(1+tanh(1.9*(ins_effect.out-0.85)));
  
  // G6P to GLU
  Real g6p_to_glu_ins_eff = 0.5*(1+tanh(1.1*(0.65-ins_effect.out))) + 0.65;
  
  // G6P to GLY
  Real g6p_to_gly_ins_eff = 2*tanh(0.55*ins_effect.out);
  
  // GLY to G6P
  Real gly_to_g6p_ins_eff = 0.985*(1+tanh(2.8*(1-ins_effect.out))+0.01);
  Real gly_to_g6p_glu_eff = 0.37*(1.15 + tanh(1.15*(1.25-glu_normed)));
  Real gly_to_g6p_gln_eff = 3.2*(tanh(0.4*(gln_glycolysis.out-0.2)));
  
  // G3P to PYR
  Real g3p_to_pyr_gln_eff = 1.18 + 1.18*(tanh(0.5*(0.7-gln_gluconeo.out)));
  
initial equation
  gln_gluconeo.out = 1;
  gln_glycolysis.out = 1;
  ins_effect.out = insulin;
  
equation

  // Hormonal effects
  gln_glycolysis.inp = glucagon;
  gln_gluconeo.inp = glucagon;
  ins_effect.inp = insulin;
  
  // Insulin release
  ins_release.G_H_inp = extra_cell_in_C[1, MPlasma.GLU];
  pancreas_insulin_inflow.k_source = ins_release.r_PIR;
  
  // Insulin clearance
  ins_clearance.k_source = -0.4 *
    (LiverConstants.input_blood_flow[1, MPlasma.INS]*extra_cell_in_C[1, MPlasma.INS] + 
     LiverConstants.input_blood_flow[2, MPlasma.INS]*extra_cell_in_C[2, MPlasma.INS] +
     ins_release.r_PIR);
     
  connect(pancreas_insulin_inflow.metabolite, extra_pool.metabolites[MPlasma.INS]);
  connect(ins_clearance.metabolite, extra_pool.metabolites[MPlasma.INS]);


  // Reactions occuring in cells and interstitium
  
  // Glycolysis reactions
  
  // GLU to G6P
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p.met_out);
  glu_to_g6p.reg_eff = 1;
  
  // GLU to G6P glucokinase
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p_glucokinase.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p_glucokinase.met_out);
  glu_to_g6p_glucokinase.reg_eff = glu_to_g6p_ins_eff;
  
  
  // G6P to G3P
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], g6p_to_g3p.met_out);
  
  // G3P to PYR
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], g3p_to_pyr.met_out);
  g3p_to_pyr.reg_eff = g3p_to_pyr_gln_eff;
  
  // PYR TO LAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], pyr_to_lac.met_out);
  
  // PYR to ACoA
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_acoa.met_inp);
  connect(locked_pool.metabolites[MLocked.ACOA], pyr_to_acoa.met_out);
  pyr_to_acoa.reg_eff = 1;
  
  
  // Glycogen cycling
  
  // G6P to GLY
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_gly.met_inp);
  connect(locked_pool.metabolites[MLocked.GLY], g6p_to_gly.met_out);
  g6p_to_gly.reg_eff = g6p_to_gly_ins_eff;
  
  // GLY to G6P
  connect(locked_pool.metabolites[MLocked.GLY], gly_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], gly_to_g6p.met_out);
  gly_to_g6p.reg_eff = gly_to_g6p_ins_eff * gly_to_g6p_glu_eff * gly_to_g6p_gln_eff;

  
  // Gluconeogenesis
  
  // G6P to GLU
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_glu.met_inp);
  connect(transfer_pool.metabolites[MChange.GLU], g6p_to_glu.met_out);
  g6p_to_glu.reg_eff = g6p_to_glu_ins_eff;

  // G3P to G6P
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], g3p_to_g6p.met_out);
  
  // GLC to G3P
  connect(transfer_pool.metabolites[MChange.GLC], glc_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], glc_to_g3p.met_out);
  
  // PYR to G3P
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], pyr_to_g3p.met_out);
  
  // LAC to PYR
  connect(transfer_pool.metabolites[MChange.LAC], lac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], lac_to_pyr.met_out);

  // AAC to PYR
  connect(transfer_pool.metabolites[MChange.AAC], aac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], aac_to_pyr.met_out);
  
  // ACoA to energy (sink)
  connect(locked_pool.metabolites[MLocked.ACOA], acoa_burn.met_inp);
   

end LiverTissue;
