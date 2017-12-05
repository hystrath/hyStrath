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
    matteoliAndMansoori

Description


SourceFiles
    matteoliAndMansooriI.H
    matteoliAndMansoori.C
    matteoliAndMansooriIO.C

\*---------------------------------------------------------------------------*/

#ifndef matteoliAndMansoori_H
#define matteoliAndMansoori_H


#include "vector.H"
#include "polyMesh.H"
#include "rdfModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class matteoliAndMansoori Declaration
\*---------------------------------------------------------------------------*/

class matteoliAndMansoori
:
	public rdfModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        word name_;

        //- temperature
        scalar T_;

        //- density
        scalar density_;

        //- no of bins
        label nBins_;

        //- maximum radial location
        scalar rMax_;

        // binwidth
        scalar binWidth_;

        //- the y axis
        scalarField g_;

        //- the x axis
        scalarField r_;

        bool reducedUnits_;

        bool readCoefficients_;

        //- coefficients
        scalar h_;
        scalar m_;
        scalar gd_;
        scalar lambda_;
        scalar alpha_;
        scalar beta_;
        scalar theta_;


    // Private Member Functions
        scalar p1
        (
            const scalar& q1,
            const scalar& q2,
            const scalar& q3
        );

        scalar p2
        (
            const scalar& q1,
            const scalar& q2,
            const scalar& q3
        );

        scalar tailPart
        (
            const scalar& y,
            const scalar& m,
            const scalar& gd,
            const scalar& lamda,
            const scalar& alpha,
            const scalar& beta
        );

        scalar initialPart
        (
            const scalar& y,
            const scalar& gd,
            const scalar& theta
        );


public:

    //- Runtime type information
    TypeName("matteoliAndMansoori");

    // Constructors

        //- Construct from components
        matteoliAndMansoori
		(
//             Time& t,
            const dictionary& dict
		);

    // Destructor

 		~matteoliAndMansoori();


    // Member Functions

        void createRdf();

        void setCoefficients();

        void setCoefficients(const dictionary& dict);


        void setRDF(radialDistribution& rdf, const Time& runTime);


//         const scalarField& g() const;
// 
//         const scalarField& r() const;
// 
//         const scalar& binWidth() const;

}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
