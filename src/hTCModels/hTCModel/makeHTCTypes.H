/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef makeHTCTypes_H
#define makeHTCTypes_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHTCTypesThermo(HTempModel, HTempType, Thermo, Table)         \
                                                                              \
    typedef HTempModel<HTempType, Thermo> HTempModel##HTempType##Thermo;          \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        HTempModel##HTempType##Thermo,                                          \
        #HTempModel"<"#HTempType","#Thermo">",                                  \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Table,                                                                \
        HTempModel##HTempType##Thermo,                                          \
        dictionary                                                            \
    );


#define makeHTCTypes(HTempModel, HTempType, Table)                       \
                                                                              \
    typedef HTempModel<HTempType> HTempModel##HTempType;                          \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        HTempModel##HTempType,                                                  \
        #HTempModel"<"#HTempType">",                                            \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Table,                                                                \
        HTempModel##HTempType,                                                  \
        dictionary                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
