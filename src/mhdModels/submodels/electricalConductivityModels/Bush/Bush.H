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
    Foam::mhd::Bush

Description
    Constant electrical conductivity coefficients for MHD model

SourceFiles
    constSigma.C

\*---------------------------------------------------------------------------*/

#ifndef Bush_H
#define Bush_H

#include "electricalConductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                            Class Bush Declaration
\*---------------------------------------------------------------------------*/

class Bush
:
    public electricalConductivityModel
{
    // Private data

        //- Temperature exponent in the Bush electrical conductivity model
        scalar n_;
        
        //- Reference temperature
        scalar T0_;
        
        //- Reference electrical conductivity
        scalar sigma0_;
        
        
public:

    //- Runtime type information
    TypeName("Bush");


    // Constructors

        //- Construct from components
        Bush(const mhdModel& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~Bush();


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
