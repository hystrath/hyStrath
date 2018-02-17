
// Write the fields needed for the re-initialisation 
// of the thermo package in runTimeEditing/hTCPropertiesDictModification.H

p.write();
U.write();
thermo.Tt().write();

forAll(Y, speciei)
{
    composition.Tv(speciei).write();
    composition.Y(speciei).write();
}

