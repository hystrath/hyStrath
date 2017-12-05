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

#include "polyForceTwoMolecules.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyForceTwoMolecules, 0);
addToRunTimeSelectionTable(polyField, polyForceTwoMolecules, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyForceTwoMolecules::setMolIds()
{

    const word molNameA(propsDict_.lookup("molSiteIdA"));
    const word molNameB(propsDict_.lookup("molSiteIdB"));

    label molIdA(findIndex(molCloud_.cP().molIds(), molNameA));
    label molIdB(findIndex(molCloud_.cP().molIds(), molNameB));

    if(molIdA == -1)
    {
        FatalErrorIn
        (
            "polyForceTwoMolecules::setSiteId()"
        )
            << "Cannot find id, molSiteIdA: " << molNameA << nl << "in idList."
            << exit(FatalError);
    }

    molSiteIdA_ = molIdA;

    if(molIdB == -1)
    {
        FatalErrorIn
        (
            "polyForceTwoMolecules::setSiteId()"
        )
            << "Cannot find id, molSiteIdB: " << molNameB << nl << "in idList."
            << exit(FatalError);
    }

    molSiteIdB_ = molIdB;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyForceTwoMolecules::polyForceTwoMolecules
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName")),
    molSiteIdA_(-1),
    molSiteIdB_(-1),
    molPointA_(propsDict_.lookup("molPointA")),
    molPointB_(propsDict_.lookup("molPointB")),
    molPositionA_(vector::zero),
    molPositionB_(vector::zero),
    cellA_(-1),
    cellB_(-1),
    molTrackingNoA_(0),
    molTrackingNoB_(0),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    rMax_(readScalar(propsDict_.lookup("rMax"))),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    binWidth_(rMax_/nBins_),
    r_(nBins_, 0.0),
    f_(nBins_, 0.0),
    mols_(nBins_, 0.0),
    timeIndex_(0),
    separation_()
{
    // read in site id (poly)
    setMolIds();

    cellA_ = mesh_.findCell(molPointA_);
    cellB_ = mesh_.findCell(molPointB_);

    forAll(r_, i)
    {
        r_[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    if(Pstream::master())
    {
    	label writeIntervalSteps = label((writeInterval_/t.deltaT().value()) + 0.5); 
    	separation_.setSize(writeIntervalSteps, 0.0);
		virA_.setSize(writeIntervalSteps, tensor::zero);
		virB_.setSize(writeIntervalSteps, tensor::zero);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyForceTwoMolecules::~polyForceTwoMolecules()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyForceTwoMolecules::createField()
{
    molTrackingNoA_ = 0;
    molPositionA_ = vector::zero;

    // pick molecule closest to starting point: molPointA_
    // (store its trackingNumber)

    if(cellA_ != -1)
    {
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellA_];

        scalar rD = GREAT;

        forAll(molsInCell, m)
        {
            polyMolecule* molI = molsInCell[m];

            scalar rIPMag = mag(molI->position() - molPointA_);

            if(rIPMag < rD)
            {
                rD = rIPMag;
                molTrackingNoA_ = molI->trackingNumber();
                molPositionA_ = molI->position();
            }
        }

        Pout << "A: tracking number: " << molTrackingNoA_ 
             << ", mol at position = " << molPositionA_ 
            << " ( reference starting point = " << molPointA_ << " ) "
            << endl;
    }

    if(Pstream::parRun())
    {
        reduce(molTrackingNoA_, sumOp<label>());
        reduce(molPositionA_, sumOp<vector>());
    }

    molTrackingNoB_ = 0;
    molPositionB_ = vector::zero;

    // pick molecule closest to starting point: molPointB_
    // (store its trackingNumber)

    if(cellB_ != -1)
    {
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellB_];

        scalar rD = GREAT;

        forAll(molsInCell, m)
        {
            polyMolecule* molI = molsInCell[m];

            scalar rIPMag = mag(molI->position() - molPointB_);

            if(rIPMag < rD)
            {
                rD = rIPMag;
                molTrackingNoB_ = molI->trackingNumber();
                molPositionB_ = molI->position();
            }
        }

        Pout << "B: tracking number: " << molTrackingNoB_ 
             << ", mol at position = " << molPositionB_ 
            << " ( reference starting point = " << molPointB_ << " ) "
            << endl;
    }

    if(Pstream::parRun())
    {
        reduce(molTrackingNoB_, sumOp<label>());
        reduce(molPositionB_, sumOp<vector>());
    }
}

void polyForceTwoMolecules::calculateField()
{
    cellA_ = -1;
    molPositionA_ = vector::zero;
    vector forceSiteA = vector::zero;
    tensor virialSiteA = tensor::zero;

    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
//             const polyMolecule::constantProperties& cP = molCloud_.constProps(mol().id());
            
            if(mol().trackingNumber() == molTrackingNoA_)
            {
                cellA_ = mol().cell();
                molPositionA_ = mol().position();
                forceSiteA = mol().a()*molCloud_.cP().mass(mol().id());
                virialSiteA = mol().rf();

                Pout<< "A: track poly mol position: " << molPositionA_ 
                    << ", current processor: " << Pstream::myProcNo()
                    << endl;
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(molPositionA_, sumOp<vector>());
        reduce(forceSiteA, sumOp<vector>());
        reduce(virialSiteA, sumOp<tensor>());
    }

    cellB_ = -1;
    molPositionB_ = vector::zero;
    vector forceSiteB = vector::zero;
    tensor virialSiteB = tensor::zero;
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(mol().trackingNumber() == molTrackingNoB_)
            {
                cellB_ = mol().cell();
                molPositionB_ = mol().position();
                forceSiteB = mol().a()*molCloud_.cP().mass(mol().id());
                virialSiteB = mol().rf();

                Pout<< "A: track poly mol position: " << molPositionB_ 
                    << ", current processor: " << Pstream::myProcNo()
                    << endl;
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(molPositionB_, sumOp<vector>());
        reduce(forceSiteB, sumOp<vector>());
        reduce(virialSiteB, sumOp<tensor>());
    }

    vector rIJ = molPositionB_ - molPositionA_;
    scalar rIJMag = mag(rIJ);

    scalar force = forceSiteB & rIJ/rIJMag;


    if(Pstream::master())
    {
        label n = label(rIJMag/binWidth_);
    
        if
        (
            (n >= 0) && (rIJMag <= rMax_)
        )
        {
            if(n == nBins_)
            {
                n--;
            }
    
            if(n < nBins_)
            {
                mols_[n] += 1.0;
                f_[n] += force;
            }
        }

        separation_[timeIndex_] = rIJMag;
        virA_[timeIndex_] = virialSiteA;
        virB_[timeIndex_] = virialSiteB;
        timeIndex_++;
    }
}

void polyForceTwoMolecules::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField f(nBins_, 0.0);

            forAll(f_, i)
            {
                if(mols_[i] > 0.0)
                {
                    f[i] = f_[i]/mols_[i];
                }
            }

            writeTimeData
            (
                timePath_,
                "polyForceTwoMolecules_"+fieldName_+"_f.xy",
                r_,
                f
            );

            const label writeIntervalSteps = ((writeInterval_/time_.deltaT().value()) + 0.5);

            scalarField timeField(writeIntervalSteps, 0.0);
    
            forAll(timeField, tT)
            {
            	timeField[tT] = time_.time().timeOutputValue()
								- writeInterval_
								+ (tT+1)*time_.deltaT().value();
            }
    
            writeTimeData
            (
                casePath_,
                "polyForceTwoMolecules"+fieldName_+"_r.xy",
                timeField,
                separation_,
                true
            );

            writeTimeData
            (
                casePath_,
                "polyForceTwoMolecules"+fieldName_+"_virA.xy",
                timeField,
                virA_,
                true
            );

            writeTimeData
            (
                casePath_,
                "polyForceTwoMolecules"+fieldName_+"_virB.xy",
                timeField,
                virB_,
                true
            );

            // reset
            timeIndex_ = 0;
            separation_ = 0.0;
            virA_=tensor::zero;
            virB_=tensor::zero;
        }
    }
}

void polyForceTwoMolecules::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyForceTwoMolecules::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyForceTwoMolecules::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
