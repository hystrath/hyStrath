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

Class
    Foam::electricalConductivityModel::constantElectricalConductivity

Description
    Constant electrical conductivity for MHD modelling

SourceFiles
    constantElectricalConductivity.C

\*---------------------------------------------------------------------------*/

#ifndef constantElectricalConductivity_H
#define constantElectricalConductivity_H

#include "electricalConductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
               Class constantElectricalConductivity Declaration
\*---------------------------------------------------------------------------*/

class constantElectricalConductivity
:
    public electricalConductivityModel
{
    // Private data

        //- Electrical conductivity scalar
        dimensionedScalar sigmaValue_;


public:

    //- Runtime type information
    TypeName("constantElectricalConductivity");


    // Constructors

        //- Construct from components
        constantElectricalConductivity
        (
            const mhdModel& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constantElectricalConductivity();


    // Member Functions

        // Edit

          //- Update electrical conductivity field
          virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
