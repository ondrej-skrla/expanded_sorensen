record GutConstants

  constant Real volume_non_cellular[MPlasma] = {11.2, 11.2, 11.2, 11.2, 9.4} "vessels + interstitium, array - in case metabolite volumes differ [dl]";
  
  constant Real volume_non_cell_exchange = 11.2 "combined interstitium and vessicular volumes for membrane transport[dl]";
  
  constant Real volume_cellular = 10.1 "[dl]";
  
  // output basal concentrations = extra cellular concentration
  constant Real basal_C_out[MPlasma] = {95.32, 11.43, 0.640, 6.03,
                                        HeartConstants.basal_C_out[MPlasma.INS]
                                        } "mg/dl";
  
  // G6P, GLY (no), G3P, PYR, ACOA
  constant Real basal_C_locked[MLocked] = {2.40, 0, 1.35, 1.35, 2.22} "[mg/dl]";
               
  // GLU, LAC, GLC (no), AAC                                         
  constant Real basal_C_transfer[MChange] = {9.00, 25.00, 0, 15.00} "[mg/dl]";
                                   
  // membrane transport constants                            
  constant Real transport_const[MChange] = {1.000, 0.545, 0, 0.624} "[dl/min]";
                                           
  // input from from arteries and from hepatic vein
  // GlU, LAC, GLC, AAC                                          
  constant Real input_blood_flow[MPlasma] = {10.1, 10.1, 10.1, 10.1, 7.2} "[dl/min]";
end GutConstants;
