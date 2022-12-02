model TestTissue
  GlucagonModule alpha_cells;
  Real ins_normed = 1;
  Real glu_normed = 2;
    
equation
  alpha_cells.ins_normed = ins_normed;
  alpha_cells.glu_normed = glu_normed;

  
end TestTissue;
