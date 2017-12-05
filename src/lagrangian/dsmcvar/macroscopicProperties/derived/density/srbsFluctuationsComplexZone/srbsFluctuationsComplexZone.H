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
    srbsFluctuationsComplexZone

Description
    Geometric field for number density

    Coarse-graining technique may be chosen (parallelised):  
        A measurement within an arbitrary cell, consists also of 
        sampling from neighbouring cells (within rCut). 

    One specie or all species handled only.

SourceFiles
    srbsFluctuationsComplexZoneI.H
    srbsFluctuationsComplexZone.C
    srbsFluctuationsComplexZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef srbsFluctuationsComplexZone_H
#define srbsFluctuationsComplexZone_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "binModel.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class srbsFluctuationsComplexZone Declaration
\*---------------------------------------------------------------------------*/

class srbsFluctuationsComplexZone
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;

//         propertyField fields_;
        autoPtr<binModel> binModel_;
        
        word fieldName_;

        labelList typeIds_;

        word regionName_;
        label regionId_;
        
        scalar resetCounter_;
        
        scalar domainLength_;

        scalarField n_;
        
        scalarField rtSin_;
        scalarField rtCos_;
        label timeIndex_;
        label nInstantSteps_;
        
        complex Rt_;
        complex RtStar_;

        List<complex> RtField_;
        List<complex> RtStarField_;
	//label counter_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("srbsFluctuationsComplexZone");

    // Constructors

        //- Construct from components
        srbsFluctuationsComplexZone
		(
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~srbsFluctuationsComplexZone();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        //- return field
//         const volScalarField& densityField() const;

        //- write field
        void writeField();

//         const propertyField& fields() const;

        void updateProperties(const dictionary& dict);
        
        void setProperties();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
