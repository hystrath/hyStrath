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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcReflectiveParticleMembranePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcReflectiveParticleMembranePatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcReflectiveParticleMembranePatch,
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcReflectiveParticleMembranePatch::readProperties()
{
    /*specularReflectionProb_ = 
        readScalar(propsDict_.lookup("reflectionProbability"));*/
}


void dsmcReflectiveParticleMembranePatch::setProperties()
{
    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
         FatalErrorIn("dsmcReflectiveParticleMembranePatch::setProperties()")
            << "Cannot have zero typeIds." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        const label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn
            (
                "dsmcReflectiveParticleMembranePatch::setProperties()"
            )
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    specularReflectionProbs_.clear();

    specularReflectionProbs_.setSize(typeIds_.size(), 0.0);

    forAll(specularReflectionProbs_, i)
    {
        specularReflectionProbs_[i] = readScalar
        (
            propsDict_.subDict("reflectionProbabilities")
                .lookup(moleculesReduced[i])
        );
    }
}


vector dsmcReflectiveParticleMembranePatch::findOriginalPosition
(
    const dsmcParcel& p, 
    const label fI
)
{
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];
    
    //Info << "patchId_ " << tab << patchId_ << tab << "neighbPatchId_ " << tab << neighbPatchId_ << endl;
    //Info << "fC " << tab << patch.faceCentres()[fB] << tab << "fCnei " << tab << patchN.faceCentres()[fB] << endl;
    
    vector orgPosition = p.position();
                
    // Limitation: the cyclic patches face centre locations must
    // only differ by one component
                
    if
    (
        patchN.faceCentres()[fI].x() != 
            patch.faceCentres()[fI].x()
     && 
        patchN.faceCentres()[fI].x() !=  
            p.tracked().currentPosition().x()
    )
    {
        orgPosition.x() = patch.faceCentres()[fI].x();
    }
    else if
    (
        patchN.faceCentres()[fI].y() != 
            patch.faceCentres()[fI].y()
      && 
        patchN.faceCentres()[fI].y() !=  
            p.tracked().currentPosition().y()      
    )
    {
        orgPosition.y() = patch.faceCentres()[fI].y();
    }
    else if
    (
        patchN.faceCentres()[fI].z() != 
            patch.faceCentres()[fI].z()
      && 
        patchN.faceCentres()[fI].z() !=  
            p.tracked().currentPosition().z() 
    )
    {
        orgPosition.z() = patch.faceCentres()[fI].z();
    }

    return orgPosition;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcReflectiveParticleMembranePatch::dsmcReflectiveParticleMembranePatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    specularReflectionProb_
    (
        0 //readScalar(propsDict_.lookup("reflectionProbability"))
    ),
    specularReflectionProbs_(),
    nReflections_(0),
    nRejections_(0)

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    
    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcReflectiveParticleMembranePatch::~dsmcReflectiveParticleMembranePatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcReflectiveParticleMembranePatch::calculateProperties()
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
            << ", ratio reflections/(reflections+rejections): " 
            << scalar(nReflections)/scalar(nReflections+nRejections)
            << endl;
    }
}



void dsmcReflectiveParticleMembranePatch::initialConfiguration()
{}


void dsmcReflectiveParticleMembranePatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    const label& iD = findIndex(typeIds_, p.typeId()); // NEW VINCENT
    
    const label faceI = p.face();

    vector nF = mesh_.faceAreas()[faceI];
    nF /= mag(nF);

    const label fA = findIndex(coupledFacesA_, faceI);
    const label fB = findIndex(coupledFacesB_, faceI);
    
    //Info << "fA " << tab << fA << tab << "fB " << tab << fB << endl;

    Random& rndGen = cloud_.rndGen();

    const scalar d = nF & p.U();

    if(d > 0) // processor boundary
    {
        //Info << "processor boundary" << endl;
        
        if(fA != -1)
        {
            const scalar pRandom = rndGen.scalar01();

            if(specularReflectionProbs_[iD] > pRandom)
            {
                //- particle specularly reflected
                const scalar Un = p.U() & nF;
    
                p.U() -= 2.0*Un*nF;

                td.switchProcessor = false;

                nReflections_++;
                
                //const vector& orgPosition = findOriginalPosition(p, fA);
                //cloud_.porousMeas().cyclicReflectionInteraction(p, orgPosition, nF);
                /*if (cloud_.measureInstantaneousMSD()) // TODO
                {
                    Info << "Incorrect calculation" << endl;
                    if (p.isTracked())
                    {
                        if (p.tracked().inPatchId() == -1)
                        {
                            const vector& orgPosition = findOriginalPosition(p, fA);
                            
                            p.tracked().updateDistanceTravelled(orgPosition);
                            
                            p.tracked().performSpecularReflectionOnDistanceTravelled(nF);
                        }
                    }
                }*/
            }
            else
            {
                //- particle passing through the membrane
                nRejections_++;
                
                const vector& orgPosition = findOriginalPosition(p, fA);
                cloud_.porousMeas().cyclicMembraneInteraction(p, orgPosition);
            }
        }
    }
    else if(d < 0) // cyclic (non-processor boundary)
    {
        if(fB != -1)
        {    
            const scalar pRandom = rndGen.scalar01();

            if(specularReflectionProbs_[iD] > pRandom)
            {
                //- particle specularly reflected
                const scalar Un = p.U() & nF;
    
                p.U() -= 2.0*Un*nF;

                nReflections_++;
                
                /*if (cloud_.measureInstantaneousMSD()) // TODO
                {
                    Info << "Incorrect calculation" << endl;
                    if (p.isTracked())
                    {
                        if (p.tracked().inPatchId() == -1)
                        {
                            const vector& orgPosition = findOriginalPosition(p, fB);
                            
                            //Info << "cur " << tab << p.tracked().currentPosition() << tab 
                            //     << "org " << tab << orgPosition << tab
                            //     << "pos " << tab << p.position() << endl;
                            
                            p.tracked().updateDistanceTravelled(orgPosition);
                            
                            //p.tracked().performSpecularReflectionOnDistanceTravelled(nF);
                        }
                    }
                }*/
            }
            else
            {
                //- particle passing through the membrane
                nRejections_++;
                
                const vector& orgPosition = findOriginalPosition(p, fB);
                cloud_.porousMeas().cyclicMembraneInteraction(p, orgPosition);
            }
        }
        /*else // NEW VINCENT
        {
            if(p.isTracked())
            {
                p.deleteTracked();
            }
        }*/
    }
}


void dsmcReflectiveParticleMembranePatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcReflectiveParticleMembranePatch::updateProperties
(
    const dictionary& newDict
)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    readProperties();
}


} // End namespace Foam

// ************************************************************************* //
