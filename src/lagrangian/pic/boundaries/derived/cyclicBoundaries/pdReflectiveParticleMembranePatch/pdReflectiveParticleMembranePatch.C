/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
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

Description

\*---------------------------------------------------------------------------*/

#include "pdReflectiveParticleMembranePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdReflectiveParticleMembranePatch, 0);

addToRunTimeSelectionTable(pdCyclicBoundary, pdReflectiveParticleMembranePatch, dictionary);


void pdReflectiveParticleMembranePatch::readProperties()
{
    p_ = (readScalar(propsDict_.lookup("reflectionProbability")));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdReflectiveParticleMembranePatch::pdReflectiveParticleMembranePatch
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    p_(readScalar(propsDict_.lookup("reflectionProbability"))),
    nReflections_(0),
    nRejections_(0)

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdReflectiveParticleMembranePatch::~pdReflectiveParticleMembranePatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void pdReflectiveParticleMembranePatch::calculateProperties()
{
    label nReflections = nReflections_;
    label nRejections = nRejections_;

    if(Pstream::parRun())
    {
        reduce(nReflections, sumOp<label>());
        reduce(nRejections, sumOp<label>());
    }

    if(nRejections > 0)
    {
        Info<< "no Reflections: " << nReflections
            << ", no Rejections: " << nRejections
            << " ratio relfections/(reflections+rejections):"
            << scalar(nReflections)/scalar(nReflections+nRejections)
            << endl;
    }
}

void pdReflectiveParticleMembranePatch::initialConfiguration()
{}



void pdReflectiveParticleMembranePatch::controlMol
(
    pdParcel& mol,
    pdParcel::trackingData& td
)
{
    const label& faceI = mol.face();

    vector nF = mesh_.faceAreas()[faceI];
//     const vector& fC  = mesh_.faceCentres()[faceI];
//     const label f = findIndex(controlPatch(), faceI);

//     bool reflect = false;

//     label f = findIndex(faces_, faceI);
    label fA = findIndex(coupledFacesA_, faceI);
    label fB = findIndex(coupledFacesB_, faceI);

    nF /= mag(nF);

    Random& rndGen = cloud_.rndGen();

    scalar d = nF & mol.U();

/*            Info<< "parcel to reflect at pos: "
                << mol.position() << ", nF: " << nF
                << " old velocity: " << mol.U()
                << " faceI: " << faceI
                << " fA: " << fA
                << " fB: " << fB
                << " fB: " << fB
                << endl; */


    if(d > 0) // processor boundary
    {
        if(fA != -1)
        {
            const scalar pRandom = rndGen.sample01<scalar>();

            if( pRandom <= p_ ) // reflect molecule
            {
                scalar Un = mol.U() & nF;

                mol.U() -= 2.0*Un*nF;

                td.switchProcessor = false;

                nReflections_++;

//                 Pout<< "Reflected!!: mol at pos: "
//                     << mol.position() << ", nF: " << nF
//                     << " tracking number: " << mol.trackingNumber()
//                     << " new velocity: " << mol.v()
//                     << endl;

            }
            else
            {

                nRejections_++;
            }
        }
    }
    else if (d < 0) // cyclic (non-processor boundary)
    {
        if(fB != -1)
        {
//             Info<< " to reflect mol at pos: "
//                 << mol.position() << ", nF: " << nF
//                 << " old velocity: " << mol.U()
//                 << " faceI: " << faceI
//                 << " fA: " << fA
//                 << " fB: " << fB
//                 << endl;

            const scalar pRandom = rndGen.sample01<scalar>();

            if( pRandom <= p_ ) // reflect molecule
            {
                scalar Un = mol.U() & nF;

                mol.U() -= 2.0*Un*nF;

                nReflections_++;

//                 Info<< "Reflected!!: mol at pos: "
//                     << mol.position()
//                     << " new velocity: " << mol.U()
//                     << endl;
            }
            else
            {

                nRejections_++;
            }
        }
    }

}



void pdReflectiveParticleMembranePatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void pdReflectiveParticleMembranePatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    readProperties();
}




} // End namespace Foam

// ************************************************************************* //
