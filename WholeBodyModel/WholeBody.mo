model WholeBody
  
  BrainTissue brain;
  MuscleTissue muscle;
  AdiposeTissue adipose;
  
  RenalTissue kidneys;
  GutTissue gut;
  LiverTissue liver;
  
  HeartTissue heart;
  OtherTissue other_tissue;
  
  GlucagonModule alpha_cells;
  
  // Glucose infusion
  VariableSource glu_infusion;
  
  // Insulin infusion
  VariableSource ins_infusion;
  
  
  // visualisation
  
  // LIVER
  Real liver_glu_to_g6p = liver.glu_to_g6p.met_inp.Q * 60;
  Real liver_glu_to_g6p_gk = liver.glu_to_g6p_glucokinase.met_inp.Q * 60;
  Real liver_g6p_to_gly = liver.g6p_to_gly.met_inp.Q * 60;
  Real liver_gly_to_g6p = liver.gly_to_g6p.met_inp.Q * 60;
  Real liver_g6p_to_g3p = liver.g6p_to_g3p.met_inp.Q * 60;
  Real liver_g3p_to_pyr = liver.g3p_to_pyr.met_inp.Q * 60;
  Real liver_pyr_to_lac = liver.pyr_to_lac.met_inp.Q * 60;
  Real liver_pyr_to_acoa = liver.pyr_to_acoa.met_inp.Q * 60;
  Real liver_acoa_sink = liver.acoa_burn.met_inp.Q * 60;
  
  // gluconeogenesis
  Real liver_lac_to_pyr = liver.lac_to_pyr.met_inp.Q * 60;
  Real liver_aac_to_pyr = liver.aac_to_pyr.met_inp.Q * 60;
  Real liver_pyr_to_g3p = liver.pyr_to_g3p.met_inp.Q * 60;
  Real liver_glc_to_g3p = liver.glc_to_g3p.met_inp.Q * 60;
  Real liver_g3p_to_g6p = liver.g3p_to_g6p.met_inp.Q * 60;
  Real liver_g6p_to_glu = liver.g6p_to_glu.met_inp.Q * 60;
  
  // membrane exchange
  Real liver_glu_out = liver.membrane_transport.extra[MChange.GLU].Q * 60;
  Real liver_lac_out = liver.membrane_transport.extra[MChange.LAC].Q * 60;
  Real liver_glc_out = liver.membrane_transport.extra[MChange.GLC].Q * 60;
  Real liver_aac_out = liver.membrane_transport.extra[MChange.AAC].Q * 60;
  
  // hormonal control
  Real liver_insulin_C = liver.extra_cell_C[MPlasma.INS];
  Real liver_glu_C = liver.extra_cell_C[MPlasma.GLU];
  Real liver_ins_release = liver.ins_release.r_PIR;
  
  
  // KIDNEYS
  Real kidneys_glu_to_g6p = kidneys.glu_to_g6p.met_inp.Q * 60;
  Real kidneys_g6p_to_g3p = kidneys.g6p_to_g3p.met_inp.Q * 60;
  Real kidneys_g3p_to_pyr = kidneys.g3p_to_pyr.met_inp.Q * 60;
  Real kidneys_pyr_to_lac = kidneys.pyr_to_lac.met_inp.Q * 60;
  Real kidneys_pyr_to_acoa = kidneys.pyr_to_acoa.met_inp.Q * 60;
  Real kidneys_acoa_sink = kidneys.acoa_burn.met_inp.Q * 60;
  Real kidneys_glu_to_lac = kidneys.glu_to_lac.met_inp.Q * 60;
  Real kidneys_uptake = kidneys.glomerulus_glu.Q * 60;
  
  
  // gluconeogenesis
  Real kidneys_lac_to_pyr = kidneys.lac_to_pyr.met_inp.Q * 60;
  Real kidneys_aac_to_pyr = kidneys.aac_to_pyr.met_inp.Q * 60;
  Real kidneys_pyr_to_g3p = kidneys.pyr_to_g3p.met_inp.Q * 60;
  Real kidneys_glc_to_g3p = kidneys.glc_to_g3p.met_inp.Q * 60;
  Real kidneys_g3p_to_g6p = kidneys.g3p_to_g6p.met_inp.Q * 60;
  Real kidneys_g6p_to_glu = kidneys.g6p_to_glu.met_inp.Q * 60;
  
  // membrane exchange
  
  Real kidneys_glu_out = -kidneys.membrane_transport.extra[MChange.GLU].Q * 60;
  Real kidneys_lac_out = kidneys.membrane_transport.extra[MChange.LAC].Q * 60;
  Real kidneys_glc_out = kidneys.membrane_transport.extra[MChange.GLC].Q * 60;
  Real kidneys_aac_out = kidneys.membrane_transport.extra[MChange.AAC].Q * 60;  
  
  Real kidneys_glu_out_mol = kidneys_glu_out / 180.16 / 70 * 1000;
  
  // MUSCLE
  Real muscle_glu_to_g6p = muscle.glu_to_g6p.met_inp.Q * 60;
  Real muscle_g6p_to_g3p = muscle.g6p_to_g3p.met_inp.Q * 60;
  Real muscle_g3p_to_pyr = muscle.g3p_to_pyr.met_inp.Q * 60;
  Real muscle_pyr_to_lac = muscle.pyr_to_lac.met_inp.Q * 60;
  Real muscle_lac_to_pyr = muscle.lac_to_pyr.met_inp.Q * 60;
  Real muscle_pyr_to_acoa = muscle.pyr_to_acoa.met_inp.Q * 60;
  Real muscle_acoa_sink = muscle.acoa_sink.met_inp.Q * 60;
  
  Real muscle_aac_release = muscle.amino_acid_source.k_source;

  Real muscle_g6p_to_gly = muscle.g6p_to_gly.met_inp.Q * 60;
  Real muscle_gly_to_g6p = muscle.gly_to_g6p.met_inp.Q * 60;
  
  Real muscle_insulin_venous = muscle.insulin_venous;
  Real muscle_insulin_interstitial = muscle.insulin;
  
  // membrane exchange
  Real muscle_glu_out = muscle.membrane_transport.extra[MChange.GLU].Q * 60;
  Real muscle_lac_out = muscle.membrane_transport.extra[MChange.LAC].Q * 60;
  Real muscle_glc_out = muscle.membrane_transport.extra[MChange.GLC].Q * 60;
  Real muscle_aac_out = muscle.membrane_transport.extra[MChange.AAC].Q * 60;
  
  
  // ADIPOSE
  Real adipose_glu_to_g6p = adipose.glu_to_g6p.met_inp.Q * 60;
  Real adipose_g6p_to_pyr = adipose.g6p_to_pyr.met_inp.Q * 60;
  Real adipose_pyr_to_lac = adipose.pyr_to_lac.met_inp.Q * 60;
  Real adipose_lac_to_pyr = adipose.lac_to_pyr.met_inp.Q * 60;
  Real adipose_pyr_to_acoa = adipose.pyr_to_acoa.met_inp.Q * 60;
  Real adipose_acoa_sink = adipose.acoa_sink.met_inp.Q * 60;
  
  Real adipose_glc_release = adipose.glycerol_cell_source.k_source;
  
  // membrane exchange
  Real adipose_glu_out = adipose.membrane_transport.extra[MChange.GLU].Q * 60;
  Real adipose_lac_out = adipose.membrane_transport.extra[MChange.LAC].Q * 60;
  Real adipose_glc_out = adipose.membrane_transport.extra[MChange.GLC].Q * 60;
  Real adipose_aac_out = adipose.membrane_transport.extra[MChange.AAC].Q * 60;
  
  Real adipose_ins_interstitial = adipose.insulin; 
  
  // BRAIN
  Real brain_glu_to_g6p = brain.glu_to_g6p.met_inp.Q * 60;
  Real brain_g6p_to_pyr = brain.g6p_to_pyr.met_inp.Q * 60;
  Real brain_pyr_to_lac = brain.pyr_to_lac.met_inp.Q * 60;
  Real brain_pyr_to_acoa = brain.pyr_to_acoa.met_inp.Q * 60;
  Real brain_acoa_sink = brain.acoa_sink.met_inp.Q * 60;
    
  // membrane exchange
  Real brain_glu_out = brain.membrane_transport.extra[MChange.GLU].Q * 60;
  Real brain_lac_out = brain.membrane_transport.extra[MChange.LAC].Q * 60;
  Real brain_glc_out = brain.membrane_transport.extra[MChange.GLC].Q * 60;
  Real brain_aac_out = brain.membrane_transport.extra[MChange.AAC].Q * 60;
  
  // GUT
  Real gut_glu_to_g6p = gut.glu_to_g6p.met_inp.Q * 60;
  Real gut_g6p_to_g3p = gut.g6p_to_g3p.met_inp.Q * 60;
  Real gut_g3p_to_pyr = gut.g3p_to_pyr.met_inp.Q * 60;
  Real gut_pyr_to_lac = gut.pyr_to_lac.met_inp.Q * 60;
  Real gut_pyr_to_acoa = gut.pyr_to_acoa.met_inp.Q * 60;
  Real gut_acoa_sink = gut.acoa_sink.met_inp.Q * 60;
  Real gut_pyr_to_aac = gut.pyr_to_aac.met_inp.Q * 60;
    
  // membrane exchange
  Real gut_glu_out = gut.membrane_transport.extra[MChange.GLU].Q * 60;
  Real gut_lac_out = gut.membrane_transport.extra[MChange.LAC].Q * 60;
  Real gut_glc_out = gut.membrane_transport.extra[MChange.GLC].Q * 60;
  Real gut_aac_out = gut.membrane_transport.extra[MChange.AAC].Q * 60;
  
  
  // HEART
  Real heart_glu_to_g6p = heart.glu_to_g6p.met_inp.Q * 60;
  Real heart_g6p_to_g3p = heart.g6p_to_g3p.met_inp.Q * 60;
  Real heart_g3p_to_pyr = heart.g3p_to_pyr.met_inp.Q * 60;
  Real heart_pyr_to_lac = heart.pyr_to_lac.met_inp.Q * 60;
  Real heart_pyr_to_acoa = heart.pyr_to_acoa.met_inp.Q * 60;
  Real heart_acoa_sink = heart.acoa_sink.met_inp.Q * 60;
    
  // membrane exchange
  Real heart_glu_out = heart.membrane_transport.extra[MChange.GLU].Q * 60;
  Real heart_lac_out = heart.membrane_transport.extra[MChange.LAC].Q * 60;
  Real heart_glc_out = heart.membrane_transport.extra[MChange.GLC].Q * 60;
  Real heart_aac_out = heart.membrane_transport.extra[MChange.AAC].Q * 60;
  
  // levels
  Real heart_glu_C = heart.extra_cell_C[MPlasma.GLU];
  Real heart_lac_C = heart.extra_cell_C[MPlasma.LAC];
  Real heart_aac_C = heart.extra_cell_C[MPlasma.AAC];
  Real heart_glc_C = heart.extra_cell_C[MPlasma.GLC];
  Real heart_ins_C = heart.extra_cell_C[MPlasma.INS];
  
  
  // OTHER TISSUE
  Real lungs_glu_in = other_tissue.lungs[MPlasma.GLU].Q * 60;
  Real lungs_lac_out = other_tissue.lungs[MPlasma.LAC].Q * 60;
  
  
  // GLUCAGON
  Real glucagon_C = alpha_cells.glucagon_C;
  
  
  // Venous blood
  Real venous_glu = (2*muscle.extra_cell_C[MPlasma.GLU] +
                     adipose.extra_cell_C[MPlasma.GLU])/3*0.84;
                     
  Real venous_ins = (2*muscle.insulin_venous +
                     adipose.insulin_venous)/3;
  
equation

  // EUGLYCEMIC HYPERINSULEMIC CLAMP
  //

  /*
  if time > 60*60 then           
    ins_infusion.k_source = 0.6*70 "mU/min";
    glu_infusion.k_source = 100*(HeartConstants.basal_C_out[MPlasma.GLU]-
                                 heart.extra_cell_C[MPlasma.GLU]);;
  else
    ins_infusion.k_source = 0;
    glu_infusion.k_source = 0;
  end if;
  */
  
  // INTRAVENOUS GLUCOSE TOLERANCE TEST
  
  if (time > 2*3600) and (time < 2*3600+3*60) then
    glu_infusion.k_source = 11666.67 "[mg/min]";
  else
    glu_infusion.k_source = 0;
  end if;
  
  ins_infusion.k_source = 0;
  
  connect(ins_infusion.metabolite, heart.arterial_pool[MPlasma.INS]);
  connect(glu_infusion.metabolite, heart.arterial_pool[MPlasma.GLU]);


  // glucagon module
  alpha_cells.insulin = heart.insulin;
  alpha_cells.glucose = heart.extra_cell_out_C[MPlasma.GLU] /
                           HeartConstants.basal_C_out[MPlasma.GLU];

  liver.glucagon = alpha_cells.gln_normed;
  
  // other tissue regulation
  other_tissue.insulin = heart.insulin;

  // connect heart to brain
  connect(heart.extra_cell_out_C, brain.extra_cell_in_C[1,:]); 

  // connect heart to muscle
  connect(heart.extra_cell_out_C, muscle.extra_cell_in_C[1,:]);
  
  // connect heart to adipose
  connect(heart.extra_cell_out_C, adipose.extra_cell_in_C[1,:]);
  
  // connect heart to kidney glomerulus
  connect(heart.arterial_pool[MPlasma.GLU], kidneys.glomerulus_glu); 
  
  //connect heart to kidneys tubules
  
  connect(heart.extra_cell_out_C, kidneys.extra_cell_in_C[1,:]);
  
  // connect heart to gut
  connect(heart.extra_cell_out_C, gut.extra_cell_in_C[1,:]);
  
  // connect heart to liver
  connect(heart.extra_cell_out_C, liver.extra_cell_in_C[1,:]);
  
  // connect gut to liver
  connect(gut.extra_cell_out_C, liver.extra_cell_in_C[2,:]);
  
  
  // connect heart to other tissue
  connect(heart.arterial_pool, other_tissue.lungs);
  connect(heart.arterial_pool, other_tissue.red_blood_cells);
  
  
  // connect brain to heart
  connect(brain.extra_cell_out_C, heart.extra_cell_in_C[1,:]);
  
  // connect muscle to heart
  connect(muscle.extra_cell_out_C, heart.extra_cell_in_C[2,:]);
  
  // connect adipose to heart
  connect(adipose.extra_cell_out_C, heart.extra_cell_in_C[3,:]);
  
  // connect kidneys to heart
  connect(kidneys.extra_cell_out_C, heart.extra_cell_in_C[4,:]);
  
  // connect liver to heart
  connect(liver.extra_cell_out_C, heart.extra_cell_in_C[5,:]);
  
end WholeBody;
