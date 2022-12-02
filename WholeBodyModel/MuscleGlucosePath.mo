model MuscleGlucosePath

  // Plasma sources and sinks - temporary
  SOURCE glucose_source(k_source=31) "[mg/min]";
  SOURCE lactate_sink(k_source=-5) "[mg/min]";

  // Pool of metabolites
  // (GLU, G6P, GLY,
  // PYR, LAC, ACOA);
  parameter Real[6] start_mass = {2380,630,300000,240,2400,200};
  CellReservoir pool(metabolite_mass(start=start_mass));
  
  // Glucose reactions
  GLUtoG6P glu_to_g6p;
  G6PtoGLY g6p_to_gly;
  GLYtoG6P gly_to_g6p;
  G6PtoPYR g6p_to_pyr;
  PYRtoLAC pyr_to_lac;
  // LACtoPYR lac_to_pyr; negligable
  PYRtoACOA pyr_to_acoa;
  ACOA_SINK acoa_sink;
  

equation
  // GLU to G6P
  connect(pool.metabolites_arr[MetEnum.GLU], glu_to_g6p.met_inp);
  connect(pool.metabolites_arr[MetEnum.G6P], glu_to_g6p.met_out);

  // G6P to GLY
  connect(pool.metabolites_arr[MetEnum.G6P], g6p_to_gly.met_inp);
  connect(pool.metabolites_arr[MetEnum.GLY], g6p_to_gly.met_out);
  
  // GLY to G6P
  connect(pool.metabolites_arr[MetEnum.GLY], gly_to_g6p.met_inp);
  connect(pool.metabolites_arr[MetEnum.G6P], gly_to_g6p.met_out);
  
  // G6P to PYR
  connect(pool.metabolites_arr[MetEnum.G6P], g6p_to_pyr.met_inp);
  connect(pool.metabolites_arr[MetEnum.PYR], g6p_to_pyr.met_out);
  
  // PYR to LAC
  connect(pool.metabolites_arr[MetEnum.PYR], pyr_to_lac.met_inp);
  connect(pool.metabolites_arr[MetEnum.LAC], pyr_to_lac.met_out);
  
  // LAC to PYR, negligable
  // connect(pool.metabolites_arr[MetEnum.LAC], lac_to_pyr.met_inp);
  // connect(pool.metabolites_arr[MetEnum.PYR], lac_to_pyr.met_out);
  
  // PYR to ACoA
  connect(pool.metabolites_arr[MetEnum.PYR], pyr_to_acoa.met_inp);
  connect(pool.metabolites_arr[MetEnum.ACOA], pyr_to_acoa.met_out);
  
  // ACoA to energy (sink)
  connect(pool.metabolites_arr[MetEnum.ACOA], acoa_sink.acoa);
  
  // constant GLU Source from plasma
  connect(pool.metabolites_arr[MetEnum.GLU], glucose_source.metabolite);
  
  // constant LAC Sink to plasma
  connect(pool.metabolites_arr[MetEnum.LAC], lactate_sink.metabolite);
  
  
end MuscleGlucosePath;