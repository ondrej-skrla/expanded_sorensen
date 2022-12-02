model RenalExcretion "Glucose excretion in kidneys"
  SingleMetabolite metabolite;
  parameter Real kidney_volume "intracellular volume [dl]";
  parameter Real thresh "extrection threshold [mg/dl]";

protected
  Real metabolite_C = metabolite.M / kidney_volume "mg/dl"; 
   
 
equation
  if metabolite_C > thresh then
    60*metabolite.Q = -330 + 0.872*metabolite_C;
  else
    60*metabolite.Q = 71 + 71 * tanh(0.011*(metabolite_C - thresh)); 
  end if;
  
end RenalExcretion;
