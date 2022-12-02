model PeakDecline "Jump to a peak and decline to a stable level"
  Real inp "indepdendent variable - peak state";
  Real out "output multiplier";
  Real decline "decline function";
  parameter Real k_dec "decline constant";
  parameter Real stable_lev "0 to 1, where to stabilize relative to peak";
  
equation
  der(decline) = ((inp - 1)*(1-stable_lev) - decline) / k_dec;
  
  out = inp - decline;

end PeakDecline;