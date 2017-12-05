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
    polyTemperatureBerendsenBinsNew

Description


SourceFiles
    polyTemperatureBerendsenBinsNewI.H
    polyTemperatureBerendsenBinsNew.C
    polyTemperatureBerendsenBinsNewIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyTemperatureBerendsenBinsNew_H
#define polyTemperatureBerendsenBinsNew_H

#include "polyStateController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyTemperatureBerendsenBinsNew Declaration
\*---------------------------------------------------------------------------*/

class polyTemperatureBerendsenBinsNew
:
	public polyStateController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- time constant
        scalar tauT_;

        // component control
        bool componentControl_;
        bool X_;
        bool Y_;
        bool Z_;

        bool peculiar_;

        autoPtr<binModel> binModel_;

        List<label> molIds_;

        scalar temperature_;
        
        bool measureFullTemperature_;


        scalar averagingCounter_;

        scalarField mass_;
        vectorField momentum_;
        vectorField velocity_;

        scalarField kE_;
        scalarField dof_;
        scalarField angularKe_;
        scalarField measuredT_;

        //- velocity rescaling term
        scalarField chi_;
        scalarField chiLin_;
        scalarField chiAng_;


        bool resetFieldsAtOutput_;
        bool output_;
        
        bool angularControl_;
        

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("polyTemperatureBerendsenBinsNew");

    // Constructors

        //- Construct from components
        polyTemperatureBerendsenBinsNew
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~polyTemperatureBerendsenBinsNew();


    // Member Functions

        //- create an initial configuration
        void initialConfiguration();

        //- control molecules at various stages of time-step
        void controlBeforeVelocityI();
        void controlBeforeMove();
        void controlBeforeForces();
        void controlDuringForces
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );        
        void controlAfterForces();
        void controlAfterVelocityII();
        
        //- calculate any required properties (end of time-step)
        void calculateProperties(); 
        
        //- output data (optional)
        void output
        (
            const fileName& fixedPathName, 
            const fileName& timePath
        );
        
        //- update properties from a modified dictionary (optional)
        void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
