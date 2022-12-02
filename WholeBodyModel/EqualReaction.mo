partial model SinkToSource "Convert one substance to another one"
  SingleMetabolite met_inp "input metabolite";
  SingleMetabolite met_out "output metabolite";
  
  parameter Real inp_out_ratio = 1 "mass ratio - from 1 gram of inp, X grams of out";
  
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

equation
  met_inp.Q = -met_out.Q / inp_out_ratio "conversion";
end SinkToSource;