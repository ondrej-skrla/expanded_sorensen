model MembraneTransport "membrane transport between extra- and intra- spaces"
  SingleMetabolite intra[MChange];
  SingleMetabolite extra[MChange];
  
  // membrane transport constant
  Real k_transport[MChange] "[dl/min]";
  parameter Real thresh[MChange] = {1000, 1000, 1000, 1000} "saturation threshold, by default very big (no threshold) mg/min";
  
  // volumes
  parameter Real intra_volume;
  parameter Real extra_volume;
  
  // Real ins_eff[MChange]; // effect of insulin on membrane transport
  
protected
  Real extra_C[MChange] = extra.M / extra_volume "concentration in the extra space [mg/dl]";
  Real intra_C[MChange] = intra.M / intra_volume "concentration in the intra space [mg/dl]";
  
equation
  
  
  60 * extra.Q = thresh .* tanh(k_transport .* (extra_C - intra_C) ./ thresh);
  intra.Q = -extra.Q;

end MembraneTransport;
