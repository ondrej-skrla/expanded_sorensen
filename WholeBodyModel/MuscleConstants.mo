record MuscleConstants

  // vessicular and interstitial volumes combined
  constant Real volume_non_cellular[MPlasma] = {70.16, 70.16, 70.16, 70.16, 6.8} "dl";
  
  constant Real volume_non_cell_exchange = 70.16 "dl";
  
  constant Real volume_cellular = 176.9 "dl";
  
  // output basal concentrations = extra cellular concentration
  constant Real basal_ins = HeartConstants.basal_C_out[MPlasma.INS]*(1-0.15);
  constant Real basal_C_out[MPlasma] = {95.62, 11.79, 0.64, 8.15, basal_ins} "mg/dl";
  
  // G6P, GLY, G3P, PYR, ACOA
  constant Real basal_C_locked[MLocked] = {2.40, 2100, 2.37, 2.37, 10.00} "[mg/dl]";
  
  // GLU, LAC, GLC (no), AAC                                         
  constant Real basal_C_transfer[MChange] = {9, 12.5, 0, 5000} "[mg/dl]";
  
  // membrane transport constants                            
  constant Real transport_const[MChange] = {0.196, 15.474, 0, 0.005} "[dl/min]";
  
  constant Real input_blood_flow[MPlasma] = {10.1, 10.1, 10.1, 10.1, 8.1} "(GLU, LAC, GLC, AAC [dl/min]";
end MuscleConstants;
