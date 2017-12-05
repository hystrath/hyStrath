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
    rdfModel

Description

SourceFiles
    rdfModelI.H
    rdfModel.C
    rdfModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef rdfModel_H
#define rdfModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Time.H"
#include "radialDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rdfModel Declaration
\*---------------------------------------------------------------------------*/

class rdfModel
{

protected:

    // Protected data

//         Time& time_;

        const dictionary& dict_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("rdfModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            rdfModel,
            dictionary,
            (
//                 Time& t,
                const dictionary& dict
            ),
            (/*t, */dict)
        );
   
    // Constructors

        //- Construct from components
        rdfModel
        (
// 			Time& t,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<rdfModel> New
        (
// 			Time& t,
            const dictionary& dict
        );


    // Destructor

        virtual ~rdfModel();


    // Member Functions

        //- create rdf
//  		virtual void createRdfFromModel() = 0;

        virtual void setRDF(radialDistribution&, const Time&) = 0;

//         virtual const scalarField& g() const = 0;
// 
//         virtual const scalarField& r() const = 0;
// 
//         virtual const label& nBins() const = 0;

//         virtual const scalar& binWidth() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
