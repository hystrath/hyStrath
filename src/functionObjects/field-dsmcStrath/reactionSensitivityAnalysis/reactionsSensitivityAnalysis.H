/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::reactionsSensitivityAnalysis

Group
    grpFieldFunctionObjects grpThermophysicalFunctionObjects

Description
    This function object creates four data files named:

    - "consumption"    :   consumption rate
    - "production"     :   destruction rate
    - "productionInt"  :   integral between dumps of the production rate
    - "consumptionInt" :   integral between dumps of the consumption rate

    The function object indicates reaction rates of creation or destruction
    of species in each reaction.


SourceFiles
    reactionsSensitivityAnalysis.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjescts_reactionsSensitivityAnalysis_H
#define functionObjescts_reactionsSensitivityAnalysis_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class reactionsSensitivityAnalysis Declaration
\*---------------------------------------------------------------------------*/

template<class chemistryType>
class reactionsSensitivityAnalysis
:
    public fvMeshFunctionObject,
    public writeFile

{
    // Private data

        //- List list for species production
        scalarListList production_;

        //- List list for species consumption
        scalarListList consumption_;

        //- List list for species production integral
        scalarListList productionInt_;

        //- List list for species consumption integral
        scalarListList consumptionInt_;

        //- Start time of integration
        scalar startTime_;

        //- End time of integration
        scalar endTime_;

        //- Word list of species
        wordList speciesNames_;

        //-Number of reactions
        label nReactions_;


         // File streams

            //- Integrated coefficients
            autoPtr<OFstream> prodFilePtr_;

            //- Moment coefficient
            autoPtr<OFstream> consFilePtr_;

            //- Drag coefficient
            autoPtr<OFstream> prodIntFilePtr_;

            //- Lift coefficient
            autoPtr<OFstream> consIntFilePtr_;



    // Private Member Functions


        //- Create file names for forces and bins
        void createFileNames();

        //- Output file header information
        void writeFileHeader(OFstream& os);

        //- Calculate production and destruction of each species
        void calculateSpeciesRR(const basicChemistryModel&);

        //- Write species production/consumption rates
        void writeSpeciesRR();


        //- Disallow default bitwise copy construct
        reactionsSensitivityAnalysis(const reactionsSensitivityAnalysis&);

        //- Disallow default bitwise assignment
        void operator=(const reactionsSensitivityAnalysis&);


public:

    //- Runtime type information
    TypeName("reactionsSensitivityAnalysis");


    // Constructors

        //- Construct from Time and dictionary
        reactionsSensitivityAnalysis
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reactionsSensitivityAnalysis();


    // Member Functions

        //- Read the reactionsSensitivityAnalysis data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the reactionsSensitivityAnalysis and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reactionsSensitivityAnalysis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
