model SlowSaturation "exponentially saturate to a stable level"
  Real inp "driving input level = saturation level";
  Real out "current level";
  parameter Real k "time constant";

initial equation
  out = inp;

equation
  der(out) = (inp - out)/k;
  

end SlowSaturation;
