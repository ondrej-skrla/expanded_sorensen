record AdiposeConstants

  constant Real volume_non_cellular[MPlasma] = {7.64, 7.64, 7.64, 7.64, 0.7} "vessels + interstitium [dl]";
  
  constant Real volume_non_cell_exchange = 7.64 "combined interstitium and vessicular volumes [dl]";
   
  constant Real volume_cellular = 19.65 "[dl]";
  
  // output basal concentrations = extra cellular concentration
  constant Real basal_ins = HeartConstants.basal_C_out[MPlasma.INS]*(1-0.15);
  constant Real basal_C_out[MPlasma] = {95.30, 11.70, 2.84, 5.48, 
                                        basal_ins} "mg/dl";
  
  // G6P, GLY, G3P, PYR, ACOA
  constant Real basal_C_locked[MLocked] = {2.40, 0, 0, 2.27, 2.22} "[mg/dl]";
  
  // GLU, LAC, GLC, AAC                                  
  constant Real basal_C_transfer[MChange] = {9.00, 25.00, 5.00, 0} "[mg/dl]";
  
  constant Real transport_const[MChange] = {0.116, 0.376, 5.093, 0} "[dl/min]";
  
  constant Real input_blood_flow[MPlasma] = {5, 5, 5, 5, 2.7} "(GLU, LAC, GLC, AAC [dl/min]";

end AdiposeConstants;
