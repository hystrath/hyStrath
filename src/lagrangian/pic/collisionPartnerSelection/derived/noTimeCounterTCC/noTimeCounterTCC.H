/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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


Description

    no-time counter with reactions

SourceFiles
    noTimeCounterTCCI.H
    noTimeCounterTCC.C
    noTimeCounterTCCIO.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeCounterTCC_H
#define noTimeCounterTCC_H


#include "collisionPartnerSelection.H"
#include "pdCloud.H"
#include "pdParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noTimeCounterTCC Declaration
\*---------------------------------------------------------------------------*/

class noTimeCounterTCC
:
    public collisionPartnerSelection
{

private:

    // Private data
        const dictionary coeffDict_;

        //- list of particles in collision groups
        List<DynamicList<pdParcel*> > collisionCellOccupancy_;

        //- list of cells in collision groups
        DynamicList<labelList> collGroupList_;

        //- desired number of particles per cell
        scalar ppC_;

        //- desired number of particles per cell
        scalar csL_;

        //- adaptation interval (number of steps)
        //label adaptInterval_;

        //- interveral counter for adaptation
        //label intervalCounter_;

        //- debug switch
        bool debug_;

        //- first run swtich
        bool firstRun_;

        //- counter for collision cell groups
        volScalarField pdColl_;

        //- visualisation of collision cell groups
        volScalarField pdRndColl_;

        //- visualisation of collision cell groups
        volScalarField pdRhoNColl_;

        scalar mfpReferenceTemperature_;

        labelList typeIds_;

        volScalarField translationalT_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField linearKE_;
        volVectorField momentum_;

        volScalarField meanFreePath_;

        List<scalarField> nParcels_;

        List<scalarField> mfp_; // mean free path

public:

    //- Runtime type information
    TypeName("noTimeCounterTCC");

    // Static data members

    // Constructors

        //- construct from components
        noTimeCounterTCC
        (
            const polyMesh& mesh,
            pdCloud& cloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~noTimeCounterTCC();

    // Member Functions

        // Access

        // Check

        // Edit
            void resetFields();

            void updateCollisionCellOccupancy();

            void initialConfiguration();

            void collide();

            void measureMeanFreePath();
        // Write





    // Member Operators

//         void operator=(const noTimeCounterTCC&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const noTimeCounterTCC&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
