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
    polyFluxController

Description

SourceFiles
    polyFluxControllerI.H
    polyFluxController.C
    polyFluxControllerIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyFluxController_H
#define polyFluxController_H

#include "IOdictionary.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"
#include "writeTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyFluxController Declaration
\*---------------------------------------------------------------------------*/

class polyFluxController
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        //- subDictionary containing the properties
        dictionary controllerDict_;

        Time& time_; 

        //- name of face zone
        word regionName_;
        label regionId_;

        scalar zoneSurfaceArea_;

        labelList internalFaces_;
        labelList processorFaces_;

        bool control_;
        bool readStateFromFile_;

        //- set all the properties below from model if required

            bool singleValueController_;
    
            // target values
            scalar density_;
            vector velocity_;
            scalar temperature_;
            scalar pressure_;
    
            tensor strainRate_;
            vector tempGradient_;
    
            bool fieldController_;
    
            //- targeted fields
            scalarField densities_;
            vectorField velocities_;
            scalarField temperatures_;
            scalarField pressures_;
    
            bool writeInTimeDir_;
            bool writeInCase_;

    // Private Member Functions

        void setFacesInfo();


public:

    //- Runtime type information
    TypeName("polyFluxController");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyFluxController,
            dictionary,
            (
                Time& t,
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (t, molCloud, dict)
        );
   
    // Constructors

        //- Construct from components
        polyFluxController
        (
			Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<polyFluxController> New
        (
			Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~polyFluxController();


    // Member Functions

        //- create an initial configuration
        virtual void initialConfiguration() = 0;


        //- calculate any required properties
        virtual void calculateProperties() = 0;


        //- control the polyMolecule from the tracking function
        virtual void controlMol
        (
            polyMolecule& mol,
            polyMolecule::trackingData& td
        ) = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;


        //- E. update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;


        void updateFluxControllerProperties(const dictionary&);

        // Access

            //- return the control zone cells
            const labelList& controlZone() const;


            label isFaceOnControlZone(const label& faceI);

            //- return the control zone name
            const word& regionName() const;

           //- return the targeted values 
            const scalar& density() const;
            scalar& density();

            const vector& velocity() const;
            vector& velocity();

            const scalar& temperature() const;
            scalar& temperature();

            const scalar& pressure() const;
            scalar& pressure();

            const tensor& strainRate() const;
            tensor& strainRate();

            const vector& tempGradient() const;
            vector& tempGradient();

            //- return the targeted fields 
            const scalarField& densityField() const;
            scalarField& densityField();

            const vectorField& velocityField() const;
            vectorField& velocityField();

            const scalarField& temperatureField() const;
            scalarField& temperatureField();

            const scalarField& pressureField() const;
            scalarField& pressureField();


            const bool& singleValueController() const;
            bool& singleValueController();

            const bool& fieldController() const;
            bool& fieldController();

            const bool& writeInTimeDir() const;
            const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
