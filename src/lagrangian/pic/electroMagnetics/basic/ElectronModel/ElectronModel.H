/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::ElectronModel

Description
    Templated PD electron fluid model class

SourceFiles
    ElectronModel.C
    ElectronModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ElectronModel_H
#define ElectronModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "pdParcel.H"

#include "fvMesh.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ElectronModel Declaration
\*---------------------------------------------------------------------------*/

class ElectronModel
{
    // Private data

protected:

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        pdCloud& cloud_;

//         //- The coefficients dictionary
//         const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("ElectronModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ElectronModel,
        dictionary,
        (
            const dictionary& dict,
            pdCloud& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        ElectronModel(pdCloud& owner);

        //- Construct from components
        ElectronModel
        (
            const dictionary& dict,
            pdCloud& owner
//             const word& type
        );


    //- Destructor
    virtual ~ElectronModel();


    //- Selector
    static autoPtr<ElectronModel> New
    (
        const dictionary& dict,
        pdCloud& owner
    );


    // Access

        //- Return the owner cloud object
//         const pdCloud& owner() const;

        //- Return non-const access to the owner cloud object
//         pdCloud& owner();

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        virtual const dictionary& coeffDict() const = 0;


    // Member Functions

        //- Flag to indicate whether model activates electron model
        virtual bool active() const = 0;

        virtual void calculatePotential() {};

        virtual void calculateRhoe() {};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "ElectronModelI.H"

#endif

// ************************************************************************* //
