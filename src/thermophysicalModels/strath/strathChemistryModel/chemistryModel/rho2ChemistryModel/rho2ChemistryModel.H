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

Class
    Foam::rho2ChemistryModel

Description
    Chemistry model for density-based thermodynamics

SourceFiles
    rho2ChemistryModelI.H
    rho2ChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef rho2ChemistryModel_H
#define rho2ChemistryModel_H

#include "basic2ChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "rho2ReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class rho2ChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class rho2ChemistryModel
:
    public basic2ChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rho2ChemistryModel(const rho2ChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const rho2ChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<rho2ReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("rho2");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rho2ChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        rho2ChemistryModel(const fvMesh& mesh);


    //- Selector
    static autoPtr<rho2ChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~rho2ChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline rho2ReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const rho2ReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rho2ChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
