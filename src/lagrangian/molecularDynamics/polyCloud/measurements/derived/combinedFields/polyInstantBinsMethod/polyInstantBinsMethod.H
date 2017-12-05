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
    polyInstantBinsMethod

Description
    The method of bins measures 1D distribution properties localised in a
    user-specified zone.

    Parallelised.

    Any binModel may be applied. So far implemented the following:

        normal, cylindrical, multi sub-zones, cuboid 

    Choice of properties to output also available.

    This is an implemented model for running transient measurements ( very costly )
    
    
SourceFiles
    polyInstantBinsMethodI.H
    polyInstantBinsMethod.C
    polyInstantBinsMethodIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyInstantBinsMethod_H
#define polyInstantBinsMethod_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyInstantBinsMethod Declaration
\*---------------------------------------------------------------------------*/

class polyInstantBinsMethod
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        fileName pathName_;
        fileName nameFile_;

        propertyField fields_;

        autoPtr<binModel> binModel_;

        label nBins_;
        
        word regionName_;
        label regionId_;
        word fieldName_;

        List<label> molIds_;


        List<DynamicList<scalar> > mols_;
        List<DynamicList<scalar> > mass_;
        List<DynamicList<vector> > mom_;
//         Field<DynamicList<vector> > velocityB_;
//         Field<DynamicList<scalar> > kE_;
//         Field<DynamicList<scalar> > angularKeSum_;
//         Field<DynamicList<scalar> > dof_;
//         Field<DynamicList<tensor> > kineticTensor_;
//         Field<DynamicList<tensor> > virialTensor_        

        // fields
//         Field<DynamicList<scalar> > N_;
        List<DynamicList<scalar> > rhoN_;
        List<DynamicList<scalar> > rhoM_;
//         Field<DynamicList<vector> > USAM_;
        List<DynamicList<vector> > UCAM_;
        List<DynamicList<scalar> > T_;
//         Field<DynamicList<tensor> > stress_;
        List<DynamicList<scalar> > p_;


        bool overideVolume_;
        List<scalar> volumes_;
        
//         scalar nAvTimeSteps_;
        
//         bool resetAtOutput_;
        
        
        
    // Private Member Functions
        
        
//             void writeToStorage();

//             bool readFromStorage();

public:

    //- Runtime type information
    TypeName("polyInstantBinsMethod");

    // Constructors

        //- Construct from components
        polyInstantBinsMethod
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyInstantBinsMethod();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        //- return field
        const propertyField& fields() const;

        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        //- write field
        void writeField();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
