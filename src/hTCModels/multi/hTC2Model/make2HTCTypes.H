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

#ifndef make2HTCTypes_H
#define make2HTCTypes_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define make2HTCTypesThermo(HTempModel, HTempType, MultiThermo, Table)        \
                                                                              \
    typedef HTempModel<HTempType, MultiThermo> HTempModel##HTempType##MultiThermo;          \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        HTempModel##HTempType##MultiThermo,                                   \
        #HTempModel"<"#HTempType","#MultiThermo">",                           \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Table,                                                                \
        HTempModel##HTempType##MultiThermo,                                   \
        dictionary                                                            \
    );


#define make2HTCTypes(HTempModel, HTempType, Table)                           \
                                                                              \
    typedef HTempModel<HTempType> HTempModel##HTempType;                      \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        HTempModel##HTempType,                                                \
        #HTempModel"<"#HTempType">",                                          \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Table,                                                                \
        HTempModel##HTempType,                                                \
        dictionary                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
