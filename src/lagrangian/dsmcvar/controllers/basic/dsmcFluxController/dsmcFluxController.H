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
    dsmcFluxController

Description
    Basic/abstract class of a flux controller


SourceFiles
    dsmcFluxControllerI.H
    dsmcFluxController.C
    dsmcFluxControllerIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcFluxController_H
#define dsmcFluxController_H

#include "IOdictionary.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "dsmcParcel.H"
#include "timeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcFluxController Declaration
\*---------------------------------------------------------------------------*/

class dsmcFluxController
{

protected:

    // Protected data

//         Time& time_;

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        Random rndGen_;

        //- subDictionary containing the properties
        dictionary controllerDict_;


        dictionary timeDict_;

        timeData time_;

        //- name of face zone
        word regionName_;
        label regionId_;
//         labelList faces_;

        scalar zoneSurfaceArea_;

        labelList internalFaces_;
        labelList processorFaces_;

        bool control_;
        bool readStateFromFile_;

        //- set all the properties below from model if required


            // target single valued controller
            bool singleValueController_;
    
            scalar density_;
            vector velocity_;
            scalar temperature_;
            scalar pressure_;
    
            tensor strainRate_;
            vector tempGradient_;

            //- targeted field controller
            bool fieldController_;

            scalarField densities_;
            vectorField velocities_;
            scalarField temperatures_;
            scalarField pressures_;

            //- writing out properties    
            bool writeInTimeDir_;
            bool writeInCase_;



    // Private Member Functions

        void setFacesInfo();


public:

    //- Runtime type information
    TypeName("dsmcFluxController");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcFluxController,
            dictionary,
            (
                Time& t,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        dsmcFluxController
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcFluxController> New
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~dsmcFluxController();


    // Member Functions

        void updateTime();

        void writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const List< Pair<scalar> >& data
        );

        //- scalar field
        void writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const scalarField& xAxis,
            const scalarField& yAxis
        );

        //- vector field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData
        );

        //- tensor field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const tensorField& yData
        );


        //- create an initial configuration
        virtual void initialConfiguration() = 0;


        //- calculate any required properties
        virtual void calculateProperties() = 0;


        //- control the dsmcParcel from the tracking function
        virtual void controlParcel
        (
            dsmcParcel& mol,
            dsmcParcel::trackingData& td
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

//             const scalar avReqDensity() const;
//             const vector avReqVelocity() const;
//             const scalar avReqTemperature() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
