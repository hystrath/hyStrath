/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    propertyField

Description

	Stores all the information for the propertyField models (state and fluxes),
    defined within the propertyFieldDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef propertyField_H
#define propertyField_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "volFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class propertyField Declaration
\*---------------------------------------------------------------------------*/

class propertyField
{

    // Private data

        const fvMesh& mesh_;
        Time& time_;
        word fieldName_;

        volScalarField sField_;
        volVectorField vField_;
        volTensorField tField_;

        scalar s_;
        vector v_;
        tensor t_;

public:

    // Constructors

        propertyField
        (
            Time& t,
            const polyMesh& mesh,
            const word& fieldName
        );



        ~propertyField();

    // Member Functions

        inline const word& fieldName() const;

        inline const volScalarField& sField() const;
        inline volScalarField& sField();
        inline const volVectorField& vField() const;
        inline volVectorField& vField();
        inline const volTensorField& tField() const;
        inline volTensorField& tField();
        inline const scalar& s() const;
        inline scalar& s();
        inline const vector& v() const;
        inline vector& v();
        inline const tensor& t() const;
        inline tensor& t() ;

    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const propertyField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "propertyFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
