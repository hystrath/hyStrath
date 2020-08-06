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

#ifndef make2ChemistrySolverTypes_H
#define make2ChemistrySolverTypes_H

#include "chemistry2Solver.H"
#include "chemistry2Model.H"

#include "no2ChemistrySolver.H"
#include "Euler2Implicit.H"
/*#include "ode2.H"*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define make2ChemistrySolverType(SS, Comp, MultiThermo)                       \
                                                                              \
    typedef SS<chemistry2Model<Comp, MultiThermo> > SS##Comp##MultiThermo;    \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##MultiThermo,                                                \
        (#SS"<" + word(Comp::typeName_())                                     \
      + "," + MultiThermo::typeName() + ">").c_str(),                         \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Comp,                                                                 \
        SS##Comp##MultiThermo,                                                \
        fvMesh                                                                \
    );


#define make2ChemistrySolverTypes(CompChemModel,MultiThermo)                  \
                                                                              \
    make2ChemistrySolverType                                                  \
    (                                                                         \
        no2ChemistrySolver,                                                   \
        CompChemModel,                                                        \
        MultiThermo                                                           \
    );                                                                        \
                                                                              \
    make2ChemistrySolverType                                                  \
    (                                                                         \
        Euler2Implicit,                                                       \
        CompChemModel,                                                        \
        MultiThermo                                                           \
    );                                                                        \
/*                                                                              \*/
/*    make2ChemistrySolverType                                                  \*/
/*    (                                                                         \*/
/*        ode2,                                                                 \*/
/*        CompChemModel,                                                        \*/
/*        MultiThermo                                                           \*/
/*    );                                                                        \*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
