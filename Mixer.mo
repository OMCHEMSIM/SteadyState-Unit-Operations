model Mixer
  import Simulator.Files.*;
  parameter Integer NOC "Number of Components", NI "Number of Input streams";
  parameter Chemsep_Database.General_Properties comp[NOC];
  parameter String outPress;
  Real outP, inP[NI];
  Real inCompMolFrac[NI, NOC] "Input stream component mol fraction", inMolFlo[NI] "Input stream Molar Flow";
  Real outCompMolFrac[NOC] "Output Stream component mol fraction", outMolFlo "Output stream Molar Flow";
  Real inTotMolEnth[NI] "Inlet molar enthalpy of each stream", outTotMolEnth "Outlet molar enthalpy";
equation
//Output Pressure
  if outPress == "Inlet_Minimum" then
    outP = min(inP);
  elseif outPress == "Inlet_Average" then
    outP = sum(inP) / NI;
  elseif outPress == "Inlet_Maximum" then
    outP = max(inP);
  end if;
//Molar Balance
  outMolFlo = sum(inMolFlo[:]);
  for i in 1:NOC loop
    outCompMolFrac[i] * outMolFlo = sum(inCompMolFrac[:, i] .* inMolFlo[:]);
  end for;
//Energy balance
  outTotMolEnth = sum(inTotMolEnth[:] .* inMolFlo[:] / outMolFlo);
end Mixer;
