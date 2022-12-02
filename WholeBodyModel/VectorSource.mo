partial model VectorSource "multi metabolite source/sink for plasma"
  parameter Integer n_inp "number of inputs";

  SingleMetabolite metabolites[MPlasma] "metabolites array";
  Real k_source[MPlasma] "speed of sources";

equation
  60 * metabolites.Q = -k_source;

end VectorSource;
