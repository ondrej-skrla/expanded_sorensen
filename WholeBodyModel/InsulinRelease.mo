model InsulinRelease "A module simulating insulin release"
  
  // variables
  Real G_H_inp "arterial glucose concentration [mg/dl]";
  Real P "potentiator [-]";
  Real P_inf "late insulin response function [-]";
  Real I "inhibitor [-]";
  Real Q "labile insulin [U]";
  Real X "peak response";
  Real S "secretion rate [pmol/min]";
  Real r_PIR "rate of pancreatic insulin release";
  
  Real S_normed;
    
  // parameters
  parameter Real S_b = 0.023 "basal insulin secretion";
  parameter Real r_b_PIR = 18.345 "TO DO, basal rate of pancreatic insulin release [mU/min]";  
  parameter Real alpha = 0.0482 "1/min";
  parameter Real beta = 0.931 "1/min";
  parameter Real K = 0.00794 "1/min";
  parameter Real Q_0 = 6.33 "";
  parameter Real gamma = 0.575 "/min";
  parameter Real beta_pir_5 = 1.11 "-";
  parameter Real M_1 = 0.00747 "1/min";
  parameter Real M_2 = 0.0958 "1/min";



initial equation
  P = P_inf;
  I = X;
  Q = (K*Q_0 + gamma*P_inf)/(K + M_1*P_inf);

equation

  S_normed = S/S_b;

  // insulin secretion
  
  S = (M_1*P_inf + M_2*(max(X-I, 0))) * Q;
  
  // peak response function
  X = InsulinPeakFunc(G_H_inp);

  //
  r_PIR = S_normed * r_b_PIR;
  
  //
  60*der(P) = alpha * (P_inf-P);

  // 
  60*der(I) = beta * (X-I);
  
  //
  60*der(Q) = K * (Q_0 - Q) + gamma*P - S;
  
  //
  P_inf = X^beta_pir_5;
  
  
end InsulinRelease;
