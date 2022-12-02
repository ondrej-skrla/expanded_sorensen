model HeartTissue

  /* LIPID CODE, NOT USED
  // Pool of venous metabolites
  SinglePool TAG_pool(metabolite_mass(start=1000000, fixed=true));
  SingleMetabolite TAG_vessel;
  Real vessel_TAG = TAG_vessel.M;
  
  // Exchange TAG block
  parameter Real k_stable = 88.41;
  TAGStabilizer tag_stabilizer(k=k_stable);
  */
  
  // specify the base model
  extends OrganModel
    (n_inp=5,
     separate=true,
     extra_cell_volume = HeartConstants.volume_non_cellular,
     intra_cell_volume = HeartConstants.volume_cellular,
                     
     input_mass_flow.in_blood_flows = HeartConstants.input_blood_flow,
     output_mass_flow.out_blood_flow = HeartConstants.output_blood_flow,
                     

     membrane_transport.extra_volume = HeartConstants.volume_non_cell_exchange,
     
     extra_cell_C(start=HeartConstants.basal_C_out, each fixed=true),
     transfer_C(start=HeartConstants.basal_C_transfer, each fixed=true),       
     locked_C(start=HeartConstants.basal_C_locked, each fixed=true)
     );
    
    
  // arterial pool, for connection of other tissues than those modeled fully as organs
  // and for kidney glomerulus
  SingleMetabolite arterial_pool[MPlasma];
  
       
  // basal masses
  parameter Real[MChange] basal_trans = HeartConstants.basal_C_transfer * intra_cell_volume "[mg]";
  parameter Real[MLocked] basal_locked = HeartConstants.basal_C_locked * intra_cell_volume "[mg]";
     
                     
  // Glycolysis reactions
  GLUtoG6P glu_to_g6p(kc_1 = 30,
                      kc_2 = 8,
                      kc_3 = 9.09,
                      basal_inp = basal_trans[MChange.GLU],
                      basal_out = basal_locked[MLocked.G6P]);
                                          
  //G6PtoG3P
  FirstOrderWithRegulation
             g6p_to_g3p(k = 3,
                        basal_inp = basal_locked[MLocked.G6P]);
  
  //G3PtoPYR
  FirstOrder
             g3p_to_pyr(k = 3,
                        basal_inp = basal_locked[MLocked.G3P]);
                        
  //LACtoPYR
  FirstOrder lac_to_pyr(k = 12,
                        basal_inp = basal_trans[MChange.LAC]);
                        
  //PYR to LAC
  FirstOrder pyr_to_lac(k = 9,
                        basal_inp = basal_locked[MLocked.PYR]);
                      
  //PYRtoACOA
  TanhInhibWithRegulation
             pyr_to_acoa(k_1 = 6,
                         k_2 = 0.5,
                         k_3 = 1,
                         basal_inp = basal_locked[MLocked.PYR],
                         basal_out = basal_locked[MLocked.ACOA]);
  
  // ACOA_SINK
  ACOA_SINK acoa_sink(k_1 = 6,
                      k_2 = 6,
                      k_3 = 0.17,
                      basal_inp = basal_locked[MLocked.ACOA]);
  
  // For visualization
  Real graph_glu_to_g6p = glu_to_g6p.met_inp.Q * 60;
  Real graph_g6p_to_g3p = g6p_to_g3p.met_inp.Q * 60;
  Real graph_g3p_to_pyr = g3p_to_pyr.met_inp.Q * 60;
  Real graph_lac_to_pyr = lac_to_pyr.met_inp.Q * 60;
  Real graph_pyr_to_acoa = pyr_to_acoa.met_inp.Q * 60;
  
  Real graph_glu_out = membrane_transport.extra[MChange.GLU].Q * 60;
  Real graph_lac_out = membrane_transport.extra[MChange.LAC].Q * 60;
  Real graph_glc_out = membrane_transport.extra[MChange.GLC].Q * 60;
  Real graph_aac_out = membrane_transport.extra[MChange.AAC].Q * 60;


  // regulation effects  
  Real insulin = extra_cell_out_C[MPlasma.INS] / 
    HeartConstants.basal_C_out[MPlasma.INS];
  
  Real glu_uptake_ins_eff = 2+2*tanh(0.2*(insulin-3.8));
    
  
initial equation
  
equation

  // connect arterial pool to extracellular pool
  connect(arterial_pool, extra_pool.metabolites);

  // Insulin effect on membrane transport
  
  membrane_transport.k_transport = HeartConstants.transport_const .*
                                   {glu_uptake_ins_eff,
                                   1, 1, 1};

  // Reactions occuring in cells and interstitium
  
  // Glycolysis reactions
  
  // GLU to G6P
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p.met_out);
  glu_to_g6p.reg_eff = glu_uptake_ins_eff;
  
  // G6P to G3P
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], g6p_to_g3p.met_out);
  g6p_to_g3p.reg_eff = glu_uptake_ins_eff;
  
  // G3P to PYR
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], g3p_to_pyr.met_out);
  
  // PYR to ACoA
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_acoa.met_inp);
  connect(locked_pool.metabolites[MLocked.ACOA], pyr_to_acoa.met_out);
  pyr_to_acoa.reg_eff = 1;
  
  // PYR to LAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], pyr_to_lac.met_out);
  
  // LAC to PYR
  connect(transfer_pool.metabolites[MChange.LAC], lac_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], lac_to_pyr.met_out);
  
  // ACoA to energy (sink)
  connect(locked_pool.metabolites[MLocked.ACOA], acoa_sink.met_inp); 
  
  

end HeartTissue;
