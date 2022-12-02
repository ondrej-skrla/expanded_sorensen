model GutTissue "includes gut, stomach, spleen, pancreas"

// specify the base model
  extends OrganModel
    (n_inp=1,
     separate=true,
     extra_cell_volume = GutConstants.volume_non_cellular,
     intra_cell_volume = GutConstants.volume_cellular,
                     
     input_mass_flow.in_blood_flows = {GutConstants.input_blood_flow},
     output_mass_flow.out_blood_flow = GutConstants.input_blood_flow,
                     

     membrane_transport.extra_volume = GutConstants.volume_non_cell_exchange,
     membrane_transport.k_transport = GutConstants.transport_const,
     membrane_transport.thresh = {20, 1000, 1000, 1000},
                     
     //hormone_mass(start=start_hormones, each fixed=true),
     //basal_hormones = RenalConstants.basal_hormones,

     extra_cell_C(start=GutConstants.basal_C_out, each fixed=true),
     transfer_C(start=GutConstants.basal_C_transfer, each fixed=true),       
     locked_C(start=GutConstants.basal_C_locked, each fixed=true)
     );
       
  // basal masses
  parameter Real[MChange] basal_trans = GutConstants.basal_C_transfer * intra_cell_volume "[mg]";
  parameter Real[MLocked] basal_locked = GutConstants.basal_C_locked * intra_cell_volume "[mg]";      
                     
  // Glycolysis reactions
  GLUtoG6P glu_to_g6p(kc_1 = 21.55,
                      kc_2 = 0.05,
                      kc_3 = 40,
                      basal_inp = basal_trans[MChange.GLU],
                      basal_out = basal_locked[MLocked.G6P]);
                                          
  //G6PtoG3P
  FirstOrder
             g6p_to_g3p(k = 20,
                        basal_inp = basal_locked[MLocked.G6P]);
  
  //G3PtoPYR
  FirstOrder
             g3p_to_pyr(k = 20,
                        basal_inp = basal_locked[MLocked.G3P]);
  //PYRtoLAC
  FirstOrder pyr_to_lac(k = 3,
                        basal_inp = basal_locked[MLocked.PYR]);
                        
  //PYRtoAAC
  FirstOrder pyr_to_aac(k = 5.6,
                        basal_inp = basal_locked[MLocked.PYR]);
                        
  //PYRtoACOA
  TanhInhibWithRegulation
            pyr_to_acoa(k_1 = 20-3-5.6,
                        k_2 = 1,
                        k_3 = 1,
                        basal_inp = basal_locked[MLocked.PYR],
                        basal_out = basal_locked[MLocked.ACOA]);
  
  // ACOA_SINK
  ACOA_SINK acoa_sink(k_1 = 20-3-5.6,
                      k_2 = 1.1,
                      k_3 = 1.5,
                      basal_inp = basal_locked[MLocked.ACOA]);
                      
                      
  // For visualization
  Real graph_glu_to_g6p = glu_to_g6p.met_inp.Q * 60;
  Real graph_g6p_to_g3p = g6p_to_g3p.met_inp.Q * 60;
  Real graph_g3p_to_pyr = g3p_to_pyr.met_inp.Q * 60;
  Real graph_pyr_to_acoa = pyr_to_acoa.met_inp.Q * 60;
  Real graph_pyr_to_lac = pyr_to_lac.met_inp.Q * 60;
  Real graph_pyr_to_aac = pyr_to_aac.met_inp.Q * 60;
  
  // membrane exchange
  
  Real graph_glu_out = membrane_transport.extra[MChange.GLU].Q * 60;
  Real graph_lac_out = membrane_transport.extra[MChange.LAC].Q * 60;
  Real graph_glc_out = membrane_transport.extra[MChange.GLC].Q * 60;
  Real graph_aac_out = membrane_transport.extra[MChange.AAC].Q * 60;
  
  
  
initial equation
  
equation

  // Reactions occuring in cells and interstitium
  
  // Glycolysis reactions
  
  // GLU to G6P
  connect(transfer_pool.metabolites[MChange.GLU], glu_to_g6p.met_inp);
  connect(locked_pool.metabolites[MLocked.G6P], glu_to_g6p.met_out);
  glu_to_g6p.reg_eff = 1;
  
  // G6P to G3P
  connect(locked_pool.metabolites[MLocked.G6P], g6p_to_g3p.met_inp);
  connect(locked_pool.metabolites[MLocked.G3P], g6p_to_g3p.met_out);
  
  // G3P to PYR
  connect(locked_pool.metabolites[MLocked.G3P], g3p_to_pyr.met_inp);
  connect(locked_pool.metabolites[MLocked.PYR], g3p_to_pyr.met_out);
  
  // PYR to LAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_lac.met_inp);
  connect(transfer_pool.metabolites[MChange.LAC], pyr_to_lac.met_out);
  
  // PYR to AAC
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_aac.met_inp);
  connect(transfer_pool.metabolites[MChange.AAC], pyr_to_aac.met_out);
  
  // PYR to ACoA
  connect(locked_pool.metabolites[MLocked.PYR], pyr_to_acoa.met_inp);
  connect(locked_pool.metabolites[MLocked.ACOA], pyr_to_acoa.met_out);
  pyr_to_acoa.reg_eff = 1;

  // ACoA to energy (sink)
  connect(locked_pool.metabolites[MLocked.ACOA], acoa_sink.met_inp); 

end GutTissue;
