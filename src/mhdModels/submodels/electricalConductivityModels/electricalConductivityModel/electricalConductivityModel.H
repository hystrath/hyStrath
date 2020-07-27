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
    Foam::electricalConductivityModel

Description


\*---------------------------------------------------------------------------*/

#ifndef electricalConductivityModel_H
#define electricalConductivityModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "mhdModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                Class electricalConductivityModel Declaration
\*---------------------------------------------------------------------------*/

class electricalConductivityModel
{

protected:

    // Protected data

        //- MHD model dictionary
        const mhdModel& dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;
        
        //- Const reference to the trans-rotational temperature field
        const volScalarField& Tt_;
        
        //- Electrical conductivity field [S/m]
        volScalarField sigma_;


public:

    //- Runtime type information
    TypeName("electricalConductivityModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            electricalConductivityModel,
            mhdModel,
            (
                const mhdModel& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        electricalConductivityModel
        (
            const mhdModel& dict,
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<electricalConductivityModel> New
    (
        const mhdModel& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~electricalConductivityModel();


    // Member Functions

        // Edit
          
          //- Update electrical conductivity field
          virtual void update() = 0;
        
        // Access

          //- Reference to the mesh
          inline const fvMesh& mesh() const
          {
              return mesh_;
          }

          //- Reference to the dictionary
          inline const mhdModel& dict() const
          {
              return dict_;
          }

          //- Const reference to the electrical conductivity field
          const volScalarField& sigma() const
          {
              return sigma_;
          }
          
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

