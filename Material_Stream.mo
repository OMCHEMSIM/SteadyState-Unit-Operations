model Material_Stream
//1-mixture, 2-Liquid, 3-Vapor
  extends Modelica.Icons.SourcesPackage;
  import Simulator.Files.*;
  //extends Thermodynamic_Packages.Raoults_Law;
  parameter Integer NOC;
  parameter Chemsep_Database.General_Properties comp[NOC];
  Real P "Pressure", T(start = 0.001) "Temperature";
  Real liqPhasMolFrac "Liquid Phase mole fraction", vapPhasMolFrac "Vapor Phase mole fraction", liqPhasMasFrac "Liquid Phase mass fraction", vapPhasMasFrac "Vapor Phase Mass fraction";
  Real totMolFlo[3] "Total molar flow", totMasFlo[3] "Total Mass Flow", MW[3](each start = 0) "Average Molecular weight of Phases";
  Real compMolFrac[3, NOC] "Component mole fraction", compMasFrac[3, NOC] "Component Mass fraction", compMolFlo[3, NOC] "Component Molar flow", compMasFlo[3, NOC] "Component Mass Fraction";
  Real phasMolSpHeat[3] "phase Molar Specific Heat", compMolSpHeat[3, NOC] "Component Molar Specific Heat";
  Real phasMolEnth[3] "Phase Molar Enthalpy", compMolEnth[3, NOC] "Component Molar Enthalpy";
  Real resMolSpHeat[3] "residual specific heat", resMolEnth[3] "residual enthalpy";
equation
//Mole Balance
  totMolFlo[1] = totMolFlo[2] + totMolFlo[3];
  compMolFrac[1, :] * totMolFlo[1] = compMolFrac[2] * totMolFlo[2] + compMolFrac[3] * totMolFlo[3];
  sum(compMolFrac[3, :]) = 1;
//sum y = 1
//component molar and mass flows
  for i in 1:NOC loop
    compMolFlo[:, i] = compMolFrac[:, i] .* totMolFlo[:];
    compMasFlo[:, i] = compMasFrac[:, i] .* totMasFlo[:];
  end for;
//phase molar and mass fractions
  liqPhasMolFrac = totMolFlo[2] / totMolFlo[1];
  vapPhasMolFrac = totMolFlo[3] / totMolFlo[1];
  liqPhasMasFrac = totMasFlo[2] / totMasFlo[1];
  vapPhasMasFrac = totMasFlo[3] / totMasFlo[1];
//Conversion between mole and mass flow
  for i in 1:NOC loop
    compMasFlo[:, i] = compMolFlo[:, i] * comp[i].MW;
  end for;
  totMasFlo[:] = totMolFlo[:] .* MW[:];
//Energy Balance
  for i in 1:NOC loop
//Specific Heat and Enthalpy calculation
    compMolSpHeat[2, i] = Thermodynamic_Functions.LiqCpId(comp[i].LiqCp, T);
    compMolSpHeat[3, i] = Thermodynamic_Functions.VapCpId(comp[i].VapCp, T);
    compMolEnth[2, i] = Thermodynamic_Functions.HLiqId(comp[i].SH, comp[i].VapCp, comp[i].HOV, comp[i].Tc, T);
    compMolEnth[3, i] = Thermodynamic_Functions.HVapId(comp[i].SH, comp[i].VapCp, comp[i].HOV, comp[i].Tc, T);
  end for;
  for i in 2:3 loop
    phasMolSpHeat[i] = sum(compMolFrac[i, :] .* compMolSpHeat[i, :]) + resMolSpHeat[i];
    phasMolEnth[i] = sum(compMolFrac[i, :] .* compMolEnth[i, :]) + resMolEnth[i];
  end for;
  phasMolSpHeat[1] = liqPhasMolFrac * phasMolSpHeat[2] + vapPhasMolFrac * phasMolSpHeat[3];
  compMolSpHeat[1, :] = compMolFrac[1, :] .* phasMolSpHeat[1];
  phasMolEnth[1] = liqPhasMolFrac * phasMolEnth[2] + vapPhasMolFrac * phasMolEnth[3];
  compMolEnth[1, :] = compMolFrac[1, :] .* phasMolEnth[1];
algorithm
  for i in 1:NOC loop
    MW[:] := MW[:] + comp[i].MW * compMolFrac[:, i];
  end for;
//for average molecular weight
end Material_Stream;
