record HeartConstants

  constant Real volume_non_cellular[MPlasma] = {14.32, 14.32, 14.32, 14.32, 9.9} "vessels + interstitium, array - in case metabolite volumes differ [dl]";
  
  constant Real volume_non_cell_exchange = 14.32 "combined interstitium and vessicular volumes for membrane transport[dl]";
  
  constant Real volume_cellular = 1.52 "[dl]";
  
  // output basal concentrations = extra cellular concentration
  constant Real basal_C_out[MPlasma] = {97.3, 10.70, 0.64, 5.48, 1.5} "mg/dl, mU/dl";
  
  // G6P, GLY (no), G3P, PYR, ACOA
  constant Real basal_C_locked[MLocked] = {2.40, 0, 0.36, 0.36, 2.22} "[mg/dl]";
               
  // GLU, LAC, GLC (no), AAC (no)                                         
  constant Real basal_C_transfer[MChange] = {9.00, 4.00, 0, 0} "[mg/dl]";
                                   
  // membrane transport constants                            
  constant Real transport_const[MChange] = {0.034, 0.448, 0, 0} "[dl/min]";
                                           
  // input from from all other organs
  // GlU, LAC, GLC, AAC                                          
  constant Real input_blood_flow[5, MPlasma] =
                {BrainConstants.input_blood_flow,
                 MuscleConstants.input_blood_flow,
                 AdiposeConstants.input_blood_flow,
                 RenalConstants.input_blood_flow,
                 LiverConstants.output_blood_flow} "[dl/min]";
                 
  constant Real[MPlasma] output_blood_flow = {43.7, 43.7, 43.7, 43.7, 31.2} "[dl/min
  ]";
  
end HeartConstants;
