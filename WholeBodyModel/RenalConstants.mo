record RenalConstants

  constant Real volume_non_cellular[MPlasma] = {6.6, 6.6, 6.6, 6.6, 5.1} "vessels + interstitium, array in case metabolite volumes differ [dl]";
  
  constant Real volume_cellular = 1.8 "[dl]";
  
  // output basal concentrations = extra cellular concentration
  
  constant Real basal_ins = HeartConstants.basal_C_out[MPlasma.INS]*(1-0.3);
  
  constant Real basal_C_out[MPlasma] = {100.17, 9.21, 0.343, 4.39,
                                        basal_ins} "mg/dl";  

  // G6P, GLY (no), G3P, PYR, ACOA
  constant Real basal_C_locked[MLocked] = {4.80, 0, 0.70, 0.40, 2.22} "[mg/dl]";
               
  // GLU, LAC, GLC, AAC                                         
  constant Real basal_C_transfer[MChange] = {144.10, 4.00, 0.05, 3.00} "[mg/dl]";
    
  constant Real volume_non_cell_exchange = 6.6 "combined interstitium and vessicular volumes for membrane transport[dl]";  
                                         
  constant Real transport_const[MChange] = {0.660, 2.876, 10.24, 7.937} "membrane transport constants [dl/min]";
                                           
  constant Real input_blood_flow[MPlasma] = {10.1, 10.1, 10.1, 10.1, 7.2} "GLU, LAC, GLC, AAC [dl/min]";
  
  constant Real excretion_thresh = 460 "when kidney excretion becomes linear [mg/dl]";
  
  //constant Real basal_hormones[Horm] = {10} "Insulin and glucagon masses at basal [?, ?]";
  
end RenalConstants;
