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
    dsmcStateController

Description
    Basic/abstract class of a state controller


SourceFiles
    dsmcStateControllerI.H
    dsmcStateController.C
    dsmcStateControllerIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcStateController_H
#define dsmcStateController_H

#include "IOdictionary.H"
// #include "Time.H"
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
                           Class dsmcStateController Declaration
\*---------------------------------------------------------------------------*/

class dsmcStateController
{

protected:

    // Protected data


        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        Random& rndGen_;

        //- subDictionary containing the properties
        dictionary controllerDict_;

        dictionary timeDict_;

        timeData time_;

        scalar timePeriod_;
        scalar initialTime_;

        //- name of control zone
        word regionName_;
        label regionId_;

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


public:

    //- Runtime type information
    TypeName("dsmcStateController");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcStateController,
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
        dsmcStateController
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcStateController> New
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~dsmcStateController();


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

        //- control dsmcParcels at different stages of the integration time-step

        virtual void controlParcelsBeforeMove() = 0;
        virtual void controlParcelsBeforeCollisions() = 0;
        virtual void controlParcelsAfterCollisions() = 0;

        //- output data
 		virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateStateControllerProperties(const dictionary&);

        // Access

            //- return the control zone cells
            const labelList& controlZone() const;

            //- return the control zone name
            const word& regionName() const;

            //- return the targeted fields 
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

            scalar avReqDensity();
            vector avReqVelocity();
            scalar avReqTemperature();
            scalar avReqPressure();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
