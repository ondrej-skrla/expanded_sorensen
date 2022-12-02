function InsulinPeakFunc "Insulin peak response function"
  
  input Real G_H "arterial glucose concentration";
  output Real X "peak response";

protected
  parameter Real beta_pir_1 = 3.27 "-";
  parameter Real beta_pir_2 = 132 "";
  parameter Real beta_pir_3 = 5.93 "-";
  parameter Real beta_pir_4 = 3.02 "-";

algorithm

  X := (G_H^beta_pir_1) / ((beta_pir_2^beta_pir_1) + (beta_pir_3*G_H^beta_pir_4));

end InsulinPeakFunc;
