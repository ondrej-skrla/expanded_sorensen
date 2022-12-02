model RenalTissue
                                                    
  // specify the base model
  extends OrganModel
    (n_inp=1,
     separate=true,
     extra_cell_volume = RenalConstants.volume_non_cellular,
     intra_cell_volume = RenalConstants.volume_cellular,
                     
     input_mass_flow.in_blood_flows = {RenalConstants.input_blood_flow},
     output_mass_flow.out_blood_flow = RenalConstants.input_blood_flow,
     membrane_transport.k_transport = RenalConstants.transport_const,           

     membrane_transport.extra_volume = RenalConstants.volume_non_cell_exchange,
     extra_cell_C(start=RenalConstants.basal_C_out, each fixed=true),
     transfer_C(start=RenalConstants.basal_C_transfer, each fixed=true),       
     locked_C(start=RenalConstants.basal_C_locked, each fixed=true)
     );
       
  // basal masses
  parameter Real[MChange] basal_trans = RenalConstants.basal_C_transfer * intra_cell_volume "[mg]";
  parameter Real[MLocked] basal_locked = RenalConstants.basal_C_locked * intra_cell_volume "[mg]";
       
  // Insulin clearance
  VariableSource ins_clearance; 
       
  // Glomerulus for glucose extraction directly from arterial blood
  SingleMetabolite glomerulus_glu;
  
  TanhInhibWithRegulation
           glu_uptake(k_1 = 20,
                      k_2 = 1,
                      k_3 = 1,
                      basal_inp = HeartConstants.basal_C_out[MPlasma.GLU] *
                                  HeartConstants.volume_non_cellular[MPlasma.GLU],
                      basal_out = basal_trans[MChange.GLU]);
       
  // Kidney renal excretion
  RenalExcretion renal_excre(kidney_volume = RenalConstants.volume_cellular,
                             thresh = RenalConstants.excretion_thresh);                   
                     
  // Glycolysis reactions
  
  // Distal tubule - shortcut directly from glucose to lactate 
  // GLU to LAC
  FirstOrderWithRegulation
             glu_to_lac(k = 4,
                        basal_inp = basal_trans[MChange.GLU]);
  
  // GLU to G6P
  GLUtoG6P   glu_to_g6p(kc_1 = 192,
                        kc_2 = 5,
                        kc_3 = 1,
                        basal_inp = basal_trans[MChange.GLU],
                        basal_out = basal_locked[MLocked.G6P]);
                               
  // G6P to G3P
  FirstOrder g6p_to_g3p(k = 16,
                        basal_inp = basal_locked[MLocked.G6P]);
                      
  // G3P to PYR
  FirstOrder g3p_to_pyr(k = 16,
                        basal_inp = basal_locked[MLocked.G3P]);
                        
  //PYRtoACOA
  TanhInhibWithRegulation
            pyr_to_acoa(k_1 = 16,
                        k_2 = 0.1,
                        k_3 = 1,
                        basal_inp = basal_locked[MLocked.PYR],
                        basal_out = basal_locked[MLocked.ACOA]);
                        
                        
  // ACOA_SINK
  ACOA_SINK acoa_burn(k_1 = 16,
                      k_2 = 3.2,
                      k_3 = 0.325,
                      basal_inp = basal_locked[MLocked.ACOA]);
  
  
  // Gluconeogenesis reactions
  // AACtoPYR
  FirstOrder aac_to_pyr(k = 11,
                        basal_inp = basal_trans[MChange.AAC]);
  
  // LACtoPYR
  FirstOrder lac_to_pyr(k = 25,
                        basal_inp = basal_trans[MChange.LAC]);
                        
  // PYRtoLAC
  FirstOrder pyr_to_lac(k = 10,
                        basal_inp = basal_locked[MLocked.PYR]);
  
  // PYRtoG3P
  FirstOrderWithRegulation
             pyr_to_g3p(k = 26,
                        basal_inp = basal_locked[MLocked.PYR]);
  
  // GLCtoG3P
  FirstOrder glc_to_g3p(k = 3,
                        basal_inp = basal_trans[MChange.GLC]);
  
  // G3PtoG6P
  FirstOrder
             g3p_to_g6p(k = 29,
                        basal_inp = basal_locked[MLocked.G3P]);
  
  // G6PtoGLU
  FirstOrderWithRegulation
             g6p_to_glu(k = 29,
                        basal_inp = basal_locked[MLocked.G6P]);

  // normed masses
  Real trans_normed[MChange] = transfer_pool.metabolites.M ./ basal_trans;
  
  // regulation effects
  Real insulin = extra_cell_C[MPlasma.INS] / 
    RenalConstants.basal_C_out[MPlasma.INS];
  
  // insulin effects on gluconeogenesis and some glucose-util reactions have a delayed action
  SlowSaturation ins_gluconeo(k=120*60);
  SlowSaturation ins_glu_util(k=120*60);
  SlowSaturation ins_glu_uptake(k=10*60);
  
  // effects on gluconeogenesis, same for G6Pase, PYR to G3P
  Real pyr_to_g3p_ins_eff = 1.07+tanh(1.4*(0.95-ins_gluconeo.out));
  
  // glucose uptake (GLUT4)
  Real glu_uptake_ins_eff = 1.1+0.2*tanh(0.7*(ins_glu_uptake.out-1.8));
  
  // glucose utilization applies to GLU to G6P, PYR to ACOA
  Real glu_util_ins_eff = 1.1+0.2*tanh(0.7*(ins_glu_util.out-1.8));
  
initial equation
  ins_gluconeo.out = insulin;
  ins_glu_util.out = insulin;
  
equation

  // Hormonal effects
  ins_gluconeo.inp = insulin;
  ins_glu_util.inp = insulin;
  ins_glu_uptake.inp = insulin;
  
  // insulin clearance
  ins_clearance.k_source =
    -0.3 * RenalConstants.input_blood_flow[MPlasma.INS]*extra_cell_in_C[1, MPlasma.INS];
  connect(ins_clearance.metabolite, extra_pool.metabolites[MPlasma.INS]);

  
  // GLUT4 translocation                                   
  // Arteries (through glomerulus) to kidney glucose inflow
  connect(glomerulus_glu, glu_uptake.met_inp);
  connect(transfer_pool.metabolites[MChange.GLU], glu_uptake.met_out);
  glu_uptake.reg_eff = glu_uptake_ins_eff;
  

  // Reactions occuring in cells and interstitium
  
  // Kidney excretion
  connect(renal_excre.metabolite, transfer_pool.metabolites[MChange.GLU]);
  
  // Glycolysis reactions

  // GLU to LAC
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], glu_to_lac.met_out);
  glu_to_lac.reg_eff = glu_util_ins_eff;
  
  // GLU to G6P
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p.met_out);
  glu_to_g6p.reg_eff = glu_util_ins_eff;
  
  // G6P to G3P
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], g6p_to_g3p.met_out);
  
  // G3P to PYR
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], g3p_to_pyr.met_out);
  
  // PYR to ACoA
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_acoa.met_inp);
  connect(locked_pool.metabolites[MLocked.ACOA], pyr_to_acoa.met_out);
  pyr_to_acoa.reg_eff = glu_util_ins_eff;
  
  // ACoA burn
  connect(locked_pool.metabolites[MLocked.ACOA], acoa_burn.met_inp);
  
  // Gluconeogenesis
  
  // G6P to GLU
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_glu.met_inp);
  connect(transfer_pool.metabolites[MChange.GLU], g6p_to_glu.met_out);
  g6p_to_glu.reg_eff = pyr_to_g3p_ins_eff;

  // G3P to G6P
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], g3p_to_g6p.met_out);
  
  // GLC to G3P
  connect(transfer_pool.metabolites[MChange.GLC], glc_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], glc_to_g3p.met_out);
  
  // PYR to G3P
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], pyr_to_g3p.met_out);
  pyr_to_g3p.reg_eff = pyr_to_g3p_ins_eff;
  
  // LAC to PYR
  connect(transfer_pool.metabolites[MChange.LAC], lac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], lac_to_pyr.met_out);
  
  // PYR to LAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], pyr_to_lac.met_out);

  // AAC to PYR
  connect(transfer_pool.metabolites[MChange.AAC], aac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], aac_to_pyr.met_out);

end RenalTissue;
