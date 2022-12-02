model OtherTissue
  // lungs
  SingleMetabolite lungs[MPlasma];
  SingleMetabolite red_blood_cells[MPlasma];
  
  
  Real insulin;  
  Real lungs_ins_eff = 2.01*tanh(0.56*insulin);

equation
  60*lungs.Q = {14 * lungs_ins_eff,
                -7 * lungs_ins_eff,
                0, 0, 0} "mg/min sink";
  60*red_blood_cells.Q = {10, -9, 0, 0, 0};
  

end OtherTissue;
