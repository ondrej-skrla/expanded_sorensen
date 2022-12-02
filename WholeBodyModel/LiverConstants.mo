record LiverConstants

  // vessels + interstitium, array - in case metabolite volumes differ 
  constant Real volume_non_cellular[MPlasma] = {13.6, 13.6, 13.6, 13.6, 11.4} "[dl]";
  
  // combined interstitium and vessicular volumes for membrane transport
  constant Real volume_non_cell_exchange = 13.6 "[dl]";
  
  constant Real volume_cellular = 11.5 "[dl]";
  
  // output basal concentrations = extra cellular concentration
  constant Real basal_ins = HeartConstants.basal_C_out[MPlasma.INS] * 1.415;
  
  constant Real basal_C_out[MPlasma] = {106.43, 9.54, 0.16, 4.57, basal_ins} "mg/dl";  
  
  // G6P, GLY, G3P, PYR, ACOA
  constant Real basal_C_locked[MLocked] = {5.00, 4000, 3.33, 3.33, 2.22} "[mg/dl]";
               
  // GLU, LAC, GLC, AAC                                         
  constant Real basal_C_transfer[MChange] = {144.1, 4.00, 0.01, 3.00} "[mg/dl]";
                                   
  // membrane transport constants                            
  constant Real transport_const[MChange] = {3.583, 4.873, 39.009, 10.826} "[dl/min]";
                                           
  // input from from arteries and from hepatic vein
  // GlU, LAC, GLC, AAC                                          
  constant Real input_blood_flow[2, MPlasma] = {{2.5, 2.5, 2.5, 2.5, 1.8},
                                                {10.1, 10.1, 10.1, 10.1, 7.2}} "[dl/min]";
                                                
  constant Real output_blood_flow[MPlasma] = LiverConstants.input_blood_flow[1] +
                                             LiverConstants.input_blood_flow[2];
  

end LiverConstants;
