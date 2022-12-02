model ACOA_SINK "Conversion of acetyl CoA to energy"
  annotation(Icon(coordinateSystem(preserveAspectRatio=true,
                  extent={{-100,-100},{100,100}}),
                graphics={Text(
                  extent={{-150,140},{150,100}},
                  textColor={0,0,255},
                  textString="%name"),
                  Rectangle(
                  extent={{-90,90},{90,-90}},
                  lineColor={255,0,0},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid)}));

  SingleMetabolite met_inp;
  
  parameter Real k_1 "basal reaction rate [mg/min]";
  parameter Real k_2 "max reaction rate, normalized [-]";
  parameter Real k_3 "steepness [-]";
  parameter Real basal_inp "basal mass of input [mg]";
  
protected
  Real inp_normed = met_inp.M/basal_inp "normalized input [-]";
  
equation
  60*met_inp.Q = k_1*k_2*tanh(k_3*inp_normed);

end ACOA_SINK;
