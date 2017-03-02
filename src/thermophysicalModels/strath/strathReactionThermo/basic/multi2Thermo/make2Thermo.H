/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

InClass
    Foam::multi2Thermo

Description
    Macros for creating basic multi2 thermo packages

\*---------------------------------------------------------------------------*/

#ifndef make2Thermo_H
#define make2Thermo_H

#include "multi2Thermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define make2ThermoTypedefs(BaseThermo,Cthermo,Mixture,Transport,Type,MultiThermo,EqnOfState,Specie)\
                                                                              \
typedef                                                                       \
    Transport                                                                 \
    <                                                                         \
        species::multiThermo                                                       \
        <                                                                     \
            MultiThermo                                                            \
            <                                                                 \
                EqnOfState                                                    \
                <                                                             \
                    Specie                                                    \
                >                                                             \
            >,                                                                \
            Type                                                              \
        >                                                                     \
    > Transport##Type##MultiThermo##EqnOfState##Specie;                            \
                                                                              \
typedef                                                                       \
    Cthermo                                                                   \
    <                                                                         \
        BaseThermo,                                                           \
        Mixture<Transport##Type##MultiThermo##EqnOfState##Specie>                  \
    > Cthermo##Mixture##Transport##Type##MultiThermo##EqnOfState##Specie;          \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##Transport##Type##MultiThermo##EqnOfState##Specie,            \
    (                                                                         \
        #Cthermo"<"#Mixture"<"                                                \
      + Transport##Type##MultiThermo##EqnOfState##Specie::typeName()               \
      + ">>"                                                                  \
    ).c_str(),                                                                \
    0                                                                         \
);


#define make2Thermo(BaseThermo,Cthermo,Mixture,Transport,Type,MultiThermo,EqnOfState,Specie)\
                                                                              \
make2ThermoTypedefs                                                            \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo,                                                                  \
    Mixture,                                                                  \
    Transport,                                                                \
    Type,                                                                     \
    MultiThermo,                                                                   \
    EqnOfState,                                                               \
    Specie                                                                    \
)                                                                             \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basic2Thermo,                                                              \
    Cthermo##Mixture##Transport##Type##MultiThermo##EqnOfState##Specie,            \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    multi2Thermo,                                                              \
    Cthermo##Mixture##Transport##Type##MultiThermo##EqnOfState##Specie,            \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Type##MultiThermo##EqnOfState##Specie,            \
    fvMesh                                                                    \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
