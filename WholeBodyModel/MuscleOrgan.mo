model MuscleTissue
  // specify the base model
  extends OrganModel
    (n_inp=1,
     separate=true,
     extra_cell_volume = MuscleConstants.volume_non_cellular,
     intra_cell_volume = MuscleConstants.volume_cellular,
                     
     input_mass_flow.in_blood_flows = {MuscleConstants.input_blood_flow},
     output_mass_flow.out_blood_flow = MuscleConstants.input_blood_flow,
                     

     membrane_transport.extra_volume = MuscleConstants.volume_non_cell_exchange,

     
     extra_cell_C(start=MuscleConstants.basal_C_out, each fixed=true),
     transfer_C(start=MuscleConstants.basal_C_transfer, each fixed=true),       
     locked_C(start=MuscleConstants.basal_C_locked, each fixed=true)
     );

     
  // basal masses
  parameter Real[MChange] basal_trans = MuscleConstants.basal_C_transfer * intra_cell_volume "[mg]";
  parameter Real[MLocked] basal_locked = MuscleConstants.basal_C_locked * intra_cell_volume "[mg]";
  
  // basal flows
  parameter Real gly_b_flow = 100;
  parameter Real glu_b_uptake = 15;
  
                     
  // Glycolysis reactions
  GLUtoG6P glu_to_g6p(kc_1 = 180,
                      kc_2 = 5,
                      kc_3 = 1,
                      basal_inp = basal_trans[MChange.GLU],
                      basal_out = basal_locked[MLocked.G6P]);
                                          
  //G6PtoG3P
  FirstOrderWithRegulation 
             g6p_to_g3p(k = 23,
                        basal_inp = basal_locked[MLocked.G6P]);
                         
  //G3PtoPYR
  FirstOrder g3p_to_pyr(k = 23,
                        basal_inp = basal_locked[MLocked.G3P]);
  
  //PYRtoACOA
  TanhInhibWithRegulation
            pyr_to_acoa(k_1 = 8,
                        k_2 = 0.1,
                        k_3 = 1,
                        basal_inp = basal_locked[MLocked.PYR],
                        basal_out = basal_locked[MLocked.ACOA]);
                        
  //PYRtoAAC
  FirstOrder pyr_to_aac(k = 4,
                        basal_inp = basal_locked[MLocked.PYR]);
                               
  //PYRtoLAC
  FirstOrder pyr_to_lac(k = 26,
                        basal_inp = basal_locked[MLocked.PYR]);
                        
  //LAC to PYR
  FirstOrder lac_to_pyr(k = 15,
                        basal_inp = basal_trans[MChange.LAC]);
                               
  //ACOA_SINK
  ACOA_SINK acoa_sink(k_1 = 8,
                      k_2 = 6.5,
                      k_3 = 0.155,
                      basal_inp = basal_locked[MLocked.ACOA]); 
  
  
  // Glycogen cycling
  // G6P to GLY
  G6PtoGLY g6p_to_gly(kc_4 = gly_b_flow,
                      kc_5 = 400000,
                      kc_6 = 0.001,
                      basal_inp = basal_locked[MLocked.G6P]);
  
  // GLY to G6P              
  MichealisMentenWithReg
           gly_to_g6p(v_max = gly_b_flow,
                      k_m = 0.01,
                      basal_inp = basal_locked[MLocked.GLY]);
                      
  // Amino acid source from proteolysis
  VariableSource amino_acid_source;
  parameter Real aac_basal_release = 25 "[mg/min]"; 
         
         
  // effects
    
  VariableSource ins_transfer "transfer from arteries to interstitium";
    
  parameter Real volume_I = 60.67 "dl";
  parameter Real k_trans_ins = 20 "min";
  parameter Real basal_ins_C = 0.568 "mU/dl";
  Real insulin_venous = extra_cell_C[MPlasma.INS] "venous insulin";
  Real insulin "interstitial insulin";
  Real insulin_normed = insulin/basal_ins_C;
    
  SlowSaturation insulin_ffa_eff(k=20*60);
  SlowSaturation insulin_oxid_eff(k=40*60);
  SlowSaturation insulin_stor_eff(k=5*60);
    
  
  // insulin effect on oxidation
  Real oxid_ins_eff = 2+1.5*tanh(0.8*(insulin_oxid_eff.out-4.4));
  Real storage_ins_eff = 25.5*tanh(0.3*insulin_stor_eff.out)-8.3;
      
  // lactate/alanine production from glucose is constant and equal to basal uptake
  
  // insulin effect on amino acid release
  Real aac_break_ins_eff = 0.9+0.6*tanh(0.4*(1.5-insulin));
  

initial equation

  insulin = basal_ins_C;

equation   

  // Insulin membrane transfer
  ins_transfer.k_source = -volume_I*(insulin_venous-insulin)/k_trans_ins;
  
  60*volume_I*der(insulin) = volume_I*(insulin_venous-insulin)/k_trans_ins -
                             3.84*insulin;                 
  
    
  insulin_ffa_eff.inp = insulin_normed;
  insulin_oxid_eff.inp = insulin_ffa_eff.out;
  insulin_stor_eff.inp = insulin_normed;

  connect(ins_transfer.metabolite, extra_pool.metabolites[MPlasma.INS]);

  // Reactions occuring in cells and interstitium_pool
  
  
  // Amino acids from proteolysis
  connect(transfer_pool.metabolites[MChange.GLC], amino_acid_source.metabolite);
  amino_acid_source.k_source = aac_basal_release * aac_break_ins_eff;
  
  // Insulin effect on membrane transport
  
  // GLUT4 translocation
  membrane_transport.k_transport = MuscleConstants.transport_const .*
                                   {oxid_ins_eff + storage_ins_eff + 1,
                                   1, 1, 1};
  
  // Glucose reactions
  
  // GLU to G6P
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p.met_out);
  glu_to_g6p.reg_eff = oxid_ins_eff + storage_ins_eff + 1;

  // G6P to GLY
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_gly.met_inp);
  connect(locked_pool.metabolites[MLocked.GLY], g6p_to_gly.met_out);
  g6p_to_gly.reg_eff = 1 + glu_b_uptake*storage_ins_eff*0.5/gly_b_flow;
  
  // GLY to G6P
  connect(locked_pool.metabolites[MLocked.GLY], gly_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], gly_to_g6p.met_out);
  gly_to_g6p.reg_eff = 1 - glu_b_uptake*storage_ins_eff*0.5/gly_b_flow;
  
   // G6P to G3P
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], g6p_to_g3p.met_out);
  g6p_to_g3p.reg_eff = (1+oxid_ins_eff)*glu_b_uptake/g6p_to_g3p.k;

  
  // G3P to PYR
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], g3p_to_pyr.met_out);
  
  // PYR to LAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], pyr_to_lac.met_out);
  
  // LAC to PYR
  connect(transfer_pool.metabolites[MChange.LAC], lac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], lac_to_pyr.met_out);
  
  // PYR to AAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_aac.met_inp);
  connect(transfer_pool.metabolites[MChange.AAC], pyr_to_aac.met_out);

  // PYR to ACoA
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_acoa.met_inp);
  connect(locked_pool.metabolites[MLocked.ACOA], pyr_to_acoa.met_out);
  pyr_to_acoa.reg_eff = oxid_ins_eff*glu_b_uptake/pyr_to_acoa.k_1;
  
  // ACoA to energy (sink)
  connect(locked_pool.metabolites[MLocked.ACOA], acoa_sink.met_inp);
  
  
  // OLD CODE
  
  // Lactate transport against the gradient
  // connect(interstitium_pool.metabolites[MetEnum.LAC], lactate_cell_sink.metabolite);
  // connect(met_vessels[MetPlasma.LAC], lactate_non_cell_source.metabolite);

  // Reaction occuring in the organ vessels
  // Endogenous TAG lipolysis, FFA goes straight to interstitium, GLC goes to output vessels
  // input is connected to arterial FFAs in the whole body module
    
  // connect(lipoprotein_pool, tg_endo_to_ffa.met_inp);
  // connect(interstitium_pool.metabolites[MetEnum.FFA], tg_endo_to_ffa.met_out_1);
  // connect(met_vessels[MetPlasma.GLC], tg_endo_to_ffa.met_out_2);
  
  
  // Fat reactions
  /*
  // FFA and glycerol to muscle TAG (TGM)
  connect(interstitium_pool.metabolites[MetEnum.FFA], ffa_to_tgm.met_inp_1);
  connect(interstitium_pool.metabolites[MetEnum.GLC], ffa_to_tgm.met_inp_2);
  connect(interstitium_pool.metabolites[MetEnum.TGM], ffa_to_tgm.met_out);
  
  // TGM to FFA and glycerol
  connect(interstitium_pool.metabolites[MetEnum.TGM], tgm_to_ffa.met_inp);
  connect(interstitium_pool.metabolites[MetEnum.FFA], tgm_to_ffa.met_out_1);
  connect(interstitium_pool.metabolites[MetEnum.GLC], tgm_to_ffa.met_out_2);
  
  // FFA to ACoA
  connect(interstitium_pool.metabolites[MetEnum.FFA], ffa_to_acoa.met_inp);
  connect(interstitium_pool.metabolites[MetEnum.ACOA], ffa_to_acoa.met_out); 
  */

  // Fat reactions
  /*                      
  FFAtoTAG ffa_to_tgm(kc_19 = 19.55,
                      kc_20 = 0.0001,
                      kc_21 = 270000,
                      basal_inp_1 = basal_masses[MetEnum.FFA],
                      basal_inp_2 = basal_masses[MetEnum.GLC]);
                      
  TAGtoFFA tgm_to_ffa(kc_19 = 20.00,
                      kc_20 = 0.001,
                      kc_21 = basal_masses[MetEnum.FFA],
                      basal = 270000);
                      
  FFAtoACOA ffa_to_acoa(kc_16 = 23.4,
                        kc_17 = 0.1,
                        kc_18 = 400,
                        basal = basal_masses[MetEnum.FFA]);
  */
  
  
  // Plasma TAG decomposition
  /*
  SingleMetabolite lipoprotein_pool "input blood lipoprotein pool [mg]";
  
  TAGtoFFA tg_endo_to_ffa(kc_19 = 8.66,
                          kc_20 = 0.01,
                          kc_21 = basal_masses[MetEnum.FFA],
                          basal = 88.41);
  */
  
  /*                        
  // Lactate membrane transport against the gradient - sink from cell, source to output vessicles
  SOURCE lactate_cell_sink(k_source = -4.19) "[mg/min]";
  SOURCE lactate_non_cell_source(k_source = 4.19) "[mg/min]";
  */


end MuscleTissue;
