
// - Output MHD fields
volVectorField F = mhd->F(U);
volScalarField Qmag = mhd->Q(U);
F.write();
Qmag.write();
//mhd->electricalConductivity().write();

