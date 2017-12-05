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

#include "polyOutputProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyOutputProperties, 0);

addToRunTimeSelectionTable(polyField, polyOutputProperties, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyOutputProperties::polyOutputProperties
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    fields_(t, mesh, "dummy"),
    accumulatedTotalLinearMomentum_(vector::zero),
    accumulatedTotalMass_(0.0),
    accumulatedTotalAngularKE_(0.0),
    accumulatedTotalLinearKE_(0.0),
    accumulatedTotalPE_(0.0),
    accumulatedTotalrDotfSum_(0.0),
    accumulatedNMols_(0.0),
    accumulatedDOFs_(0.0),
    averageTemperature_(0.0),
    averagePressure_(0.0),
    nAvTimeSteps_(0.0)
{
    const scalarField& cellVols = mesh_.cellVolumes();
    
    meshVolume_ = sum(cellVols);
    
    if (Pstream::parRun())
    {
        reduce(meshVolume_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyOutputProperties::~polyOutputProperties()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyOutputProperties::createField()
{}


void polyOutputProperties::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    vector singleStepTotalLinearMomentum(vector::zero);
    vector singleStepTotalAngularMomentum(vector::zero);
    scalar singleStepMaxVelocityMag = 0.0;
    scalar singleStepTotalMass = 0.0;
    scalar singleStepTotalLinearKE = 0.0;
    scalar singleStepTotalAngularKE = 0.0;
    scalar singleStepTotalPE = 0.0;
    scalar singleStepTotalrDotf = 0.0;

    label singleStepNMols = molCloud_.size();
    label singleStepDOFs = 0;

    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            label molId = mol().id();
            scalar molMass(molCloud_.cP().mass(molId));
            singleStepTotalMass += molMass;
        }
    
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            label molId = mol().id();
    
            scalar molMass(molCloud_.cP().mass(molId));

            const vector& molV(mol().v());

            vector molPiGlobal = vector::zero;
            vector molOmega = vector::zero;
            const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molId));

            if(!molCloud_.cP().pointMolecule(molId))
            {
                molOmega = inv(molMoI) & mol().pi();
                molPiGlobal = mol().Q() & mol().pi();
            }

            singleStepTotalLinearMomentum += molV * molMass;
            singleStepTotalAngularMomentum += molPiGlobal;

            if(mag(molV) > singleStepMaxVelocityMag)
            {
                singleStepMaxVelocityMag = mag(molV);
            }
    
            singleStepTotalLinearKE += 0.5*molMass*magSqr(molV);
            singleStepTotalAngularKE += 0.5*(molOmega & molMoI & molOmega);
            singleStepTotalPE += mol().potentialEnergy();
            singleStepTotalrDotf += tr(mol().rf());
            singleStepDOFs += molCloud_.cP().degreesOfFreedom(molId);
        }
    }

    if (Pstream::parRun())
    {
        reduce(singleStepTotalLinearMomentum, sumOp<vector>());
        reduce(singleStepTotalAngularMomentum, sumOp<vector>());
        reduce(singleStepMaxVelocityMag, maxOp<scalar>());
        reduce(singleStepTotalMass, sumOp<scalar>());
        reduce(singleStepTotalLinearKE, sumOp<scalar>());
        reduce(singleStepTotalAngularKE, sumOp<scalar>());
        reduce(singleStepTotalPE, sumOp<scalar>());
        reduce(singleStepTotalrDotf, sumOp<scalar>());
        reduce(singleStepNMols, sumOp<label>());
        reduce(singleStepDOFs, sumOp<label>());
    }

    if (singleStepNMols)
    {
        Info<< "Number of molecules in system = "
            << singleStepNMols << nl
            << "Overall number density = "
            << singleStepNMols/meshVolume_ << nl
            << "Overall mass density = "
            << singleStepTotalMass/meshVolume_ << nl
            << "Average linear momentum per molecule = "
            << singleStepTotalLinearMomentum/singleStepNMols << ' '
            << mag(singleStepTotalLinearMomentum)/singleStepNMols << nl
            << "Average angular momentum per molecule = "
            << singleStepTotalAngularMomentum << ' '
            << mag(singleStepTotalAngularMomentum)/singleStepNMols << nl
            << "Maximum |velocity| = "
            << singleStepMaxVelocityMag << nl
            << "Average linear KE per molecule = "
            << singleStepTotalLinearKE/singleStepNMols << nl
            << "Average angular KE per molecule = "
            << singleStepTotalAngularKE/singleStepNMols << nl
            << "Average PE per molecule = "
            << singleStepTotalPE/singleStepNMols << nl
            << "Average TE per molecule = "
            <<
            (
                singleStepTotalLinearKE
            + singleStepTotalAngularKE
            + singleStepTotalPE
            )
            /singleStepNMols
            << endl;
    }
    else
    {
        Info<< "No molecules in system" << endl;
    }

    accumulatedTotalLinearMomentum_ += singleStepTotalLinearMomentum;
    accumulatedTotalMass_ += singleStepTotalMass;
    accumulatedTotalLinearKE_ += singleStepTotalLinearKE;
    accumulatedTotalAngularKE_ += singleStepTotalAngularKE;
    accumulatedTotalPE_ += singleStepTotalPE;
    accumulatedTotalrDotfSum_ += singleStepTotalrDotf;
    accumulatedNMols_ += singleStepNMols;
    accumulatedDOFs_ += singleStepDOFs;

}

void polyOutputProperties::writeField()
{
    const Time& runTime = time_.time();

    if (runTime.outputTime())
    {
    	const scalar& nAvTimeSteps = nAvTimeSteps_; 

        if (accumulatedNMols_)
        {
            Info << "calculating averages" << endl;
    
            averageTemperature_ =
            (
                2.0/(molCloud_.redUnits().kB() * accumulatedDOFs_)
                *
                (
                    accumulatedTotalLinearKE_ + accumulatedTotalAngularKE_
                    -
                    0.5*magSqr(accumulatedTotalLinearMomentum_)/accumulatedTotalMass_
                )
            );
    
            averagePressure_ =
            (
                (
                    (accumulatedNMols_/nAvTimeSteps)
                    *
                    molCloud_.redUnits().kB() * averageTemperature_
                    +
                    accumulatedTotalrDotfSum_/(6.0 * nAvTimeSteps)
                )
                /
                meshVolume_
            );
    
            Info << "----------------------------------------" << nl
                << "Averaged properties" << nl
                << "Average |velocity| = "
                << mag(accumulatedTotalLinearMomentum_)/accumulatedTotalMass_ << nl
                << "Average temperature = " << averageTemperature_ << nl
                << "Average pressure = " << averagePressure_ << nl
                << "----------------------------------------" << endl;
        }
        else
        {
            Info<< "Not averaging temperature and pressure: "
                << "no molecules in system" << endl;
        }

        //-reset
        accumulatedTotalLinearMomentum_ = vector::zero;
        accumulatedTotalMass_ = 0.0;
        accumulatedTotalLinearKE_ = 0.0;
        accumulatedTotalAngularKE_ = 0.0;
        accumulatedTotalPE_ = 0.0;
        accumulatedTotalrDotfSum_ = 0.0;
        accumulatedNMols_ = 0.0;
        accumulatedDOFs_ = 0.0;
        nAvTimeSteps_ = 0.0;
    }
}

void polyOutputProperties::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyOutputProperties::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyOutputProperties::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
