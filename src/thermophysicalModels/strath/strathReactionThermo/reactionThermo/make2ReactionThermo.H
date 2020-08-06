/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef make2ReactionThermo_H
#define make2ReactionThermo_H

#include "addToRunTimeSelectionTable.H"
#include "Specie2Mixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define make2ReactionMixtureThermo(BaseThermo,CThermo,MixtureThermo,Mixture,ThermoPhys)   \
                                                                              \
typedef MixtureThermo                                                         \
<                                                                             \
    CThermo,                                                                  \
    Specie2Mixture                                                            \
    <                                                                         \
        Mixture                                                               \
        <                                                                     \
            ThermoPhys                                                        \
        >                                                                     \
    >                                                                         \
> MixtureThermo##Mixture##ThermoPhys;                                         \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    (#MixtureThermo"<"#Mixture"<" + ThermoPhys::typeName() + ">>").c_str(),   \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basic2Thermo,                                                             \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    multi2Thermo,                                                             \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    CThermo,                                                                  \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    fvMesh                                                                    \
);


#define make2ReactionThermo(BaseThermo,CThermo,MixtureThermo,Mixture,Transport,Type,MultiThermo,EqnOfState,Specie) \
                                                                              \
typedef                                                                       \
    Transport                                                                 \
    <                                                                         \
        species::multiThermo                                                  \
        <                                                                     \
            MultiThermo                                                       \
            <                                                                 \
                EqnOfState                                                    \
                <                                                             \
                    Specie                                                    \
                >                                                             \
            >,                                                                \
            Type                                                              \
        >                                                                     \
    > Transport##Type##MultiThermo##EqnOfState##Specie;                       \
                                                                              \
make2ReactionMixtureThermo                                                    \
(                                                                             \
    BaseThermo,                                                               \
    CThermo,                                                                  \
    MixtureThermo,                                                            \
    Mixture,                                                                  \
    Transport##Type##MultiThermo##EqnOfState##Specie                          \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
