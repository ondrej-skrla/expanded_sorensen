model OrganModel

  // parameters and inputs
  parameter Integer n_inp "number of input organs";
  parameter Real intra_cell_volume "dl";
  parameter Real extra_cell_volume[MPlasma] "[dl]";

  // metabolite pools
  ExtraPool extra_pool "extracellular metabolites";
  TransferPool transfer_pool "intracellular metabolites";
  LockedPool locked_pool "intracellular metabolites";
  
  Real extra_cell_C[MPlasma] = extra_pool.metabolites.M ./ extra_cell_volume;
  Real transfer_C[MChange] = transfer_pool.metabolites.M / intra_cell_volume;
  Real locked_C[MLocked] = locked_pool.metabolites.M / intra_cell_volume;
  
  // inflow and outflow of metabolites via blood
  RealInput extra_cell_in_C[n_inp, MPlasma] "concentrations from input organs [mg/dl]";
  RealOutput extra_cell_out_C[MPlasma] = extra_pool.metabolites.M ./ extra_cell_volume "output concentration [mg/dl]";
    
  InputMassFlow input_mass_flow(n_inp=n_inp, vessel_in_C = extra_cell_in_C);
  OutputMassFlow output_mass_flow(n_inp=n_inp, vessel_out_C = extra_cell_out_C);
  
  
  
  // transport from venous to interstitium
  parameter Boolean separate = true "true if vessels and interstitium separated";
  MembraneTransport membrane_transport(intra_volume = intra_cell_volume);
  
  // connectors for exchange between intracellular and extracellular
  SingleMetabolite exchange_extra_cell[MChange];
  
initial equation
  //der(extra_cell_out_C) = fill(0, size(extra_cell_out_C, 1));  
  
equation

  // vessel exchange
  connect(extra_pool.metabolites, input_mass_flow.metabolites);
  connect(extra_pool.metabolites, output_mass_flow.metabolites);
  
  // connect selected metabolites from plasma to membrane flows  
  connect(extra_pool.metabolites[MPlasma.GLU], exchange_extra_cell[MChange.GLU]);
  connect(extra_pool.metabolites[MPlasma.LAC], exchange_extra_cell[MChange.LAC]);
  connect(extra_pool.metabolites[MPlasma.GLC], exchange_extra_cell[MChange.GLC]);
  connect(extra_pool.metabolites[MPlasma.AAC], exchange_extra_cell[MChange.AAC]);
  
  // if separate, model membrane transport, if not, join directly
  if separate then
     connect(transfer_pool.metabolites, membrane_transport.intra);
     connect(exchange_extra_cell, membrane_transport.extra);
  else
     connect(exchange_extra_cell, transfer_pool);
  end if;

end OrganModel;
