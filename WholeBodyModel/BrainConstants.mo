record BrainConstants
  
  // vessels + interstitium, array in case metabolite volumes differ   
  constant Real volume_non_cellular[MPlasma] = {8, 8, 8, 8, 2.5} "[dl]";
  
  constant Real volume_non_cell_exchange = 8 "[dl]";  
  
  constant Real volume_cellular = 8.6 "[dl]";
  
  // output basal concentrations = extra cellular concentration
  constant Real basal_C_out[MPlasma] = {81.84, 11.65, 0.64, 5.48,
                                        HeartConstants.basal_C_out[MPlasma.INS]
                                        } "mg/dl";
  
  // G6P, GLY (no), G3P (no), PYR, ACOA
  constant Real basal_C_locked[MLocked] = {4.80, 0, 0, 1.09, 3.00} "[mg/dl]";
  
  // GLU, LAC, GLC (no), AAC (no)                                  
  constant Real basal_C_transfer[MChange] = {18.00, 25.00, 0, 0} "[mg/dl]";
                                         
  constant Real transport_const[MChange] = {5.000, 0.419, 0, 0} "membrane transport constants [dl/min]";
  
  constant Real input_blood_flow[MPlasma] = {5.9, 5.9, 5.9, 5.9, 4.5} "(GLU, LAC, GLC, AAC [dl/min]";
  
end BrainConstants;
