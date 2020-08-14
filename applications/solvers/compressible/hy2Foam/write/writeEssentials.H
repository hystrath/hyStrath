
// Write the fields needed for the re-initialisation
// of the thermo package in runTimeEditing/hTCPropertiesDictModification.H

p.write();
U.write();
thermo.T().write();

forAll(Y, speciei)
{
    if (not (downgradeSingleT or downgradeSingleTv))
    {
        composition.Tv(speciei).write();
    }
    composition.Y(speciei).write();
}

if (downgradeSingleTv)
{
    thermo.Tv().write();
}

