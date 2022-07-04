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

Description

\*---------------------------------------------------------------------------*/

#include "polyPropertyFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyPropertyFields, 0);

addToRunTimeSelectionTable(polyField, polyPropertyFields, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPropertyFields::polyPropertyFields
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
    N_
    (
        IOobject
        (
            fieldName_+"_N",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    rhoN_
    (
        IOobject
        (
            fieldName_+"_rhoN",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            fieldName_+"_rhoM",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    ),
    USAM_
    (
        IOobject
        (
            fieldName_+"_U_SAM",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    UCAM_
    (
        IOobject
        (
            fieldName_+"_U_CAM",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    T_
    (
        IOobject
        (
            fieldName_+"_T",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    p_
    (
        IOobject
        (
            fieldName_+"_p",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    stress_
    (
        IOobject
        (
            fieldName_+"_stress",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("0.0", dimless, tensor::zero)
    ),
//     volumes_(mesh_.nCells(), scalar(0.0)),
    mols_(mesh_.nCells(), 0.0),
    mass_(mesh_.nCells(), 0.0),
    mom_(mesh_.nCells(), vector::zero),
    velocityB_(mesh_.nCells(), vector::zero),
    kE_(mesh_.nCells(), 0.0),
    dof_(mesh_.nCells(), 0.0),
    angularKe_(mesh_.nCells(), 0.0),
    kineticTensor_(mesh_.nCells(), tensor::zero),
    virialTensor_(mesh_.nCells(), tensor::zero),

    molsV_(mesh_.nCells(), 0.0),
    massV_(mesh_.nCells(), 0.0),
    momV_(mesh_.nCells(), vector::zero),
    velocity_(mesh_.nCells(), vector::zero),
    angularSpeed_(mesh_.nCells(), vector::zero),
    angularVelocity_(mesh_.nCells(), vector::zero),

    molIds_(),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
//     coarseGrainedMeasurements_(false)
{

    // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));

//     if (propsDict_.found("coarseGrainedMeasurements"))
//     {
//         coarseGrainedMeasurements_ = Switch(propsDict_.lookup("coarseGrainedMeasurements"));
//     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPropertyFields::~polyPropertyFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial condition
void polyPropertyFields::createField()
{
//     setVolumeMeasurements();

    Info << "Initialising polyPropertyFields field" << endl;

    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const scalarField& volField = mesh_.cellVolumes();

    const scalar& kB = molCloud_.redUnits().kB();

    forAll(cellOccupancy, cell)
    {
        const List<polyMolecule*>& molsInCell = cellOccupancy[cell];

        vector velocity(vector::zero);
        vector angularVelocity(vector::zero);

        {
            scalar mols = 0.0;
            scalar mass = 0.0;
            vector mom = vector::zero;
            vector angularSpeed = vector::zero;

            forAll(molsInCell, mIC)
            {
                polyMolecule* molI = molsInCell[mIC];

                if(findIndex(molIds_, molI->id()) != -1)
                {
                    mols += 1.0;

//                     const polyMolecule::constantProperties& constProp = molCloud_.constProps(molI->id());
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    mass += massI;
                    mom += massI*molI->v();

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularSpeed += molOmega;
                }
            }

            if(mass > 0.0)
            {
                velocity = mom/mass;
                angularVelocity = angularSpeed/mols;
            }
        }

        scalar mols = 0.0;
        scalar mass = 0.0;
        vector momentum = vector::zero;
        scalar kE = 0.0;
        label dof = 0;
        scalar angularKeSum = 0.0;

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            if(findIndex(molIds_, molI->id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(molI->id());

                mass += massI;
                mols += 1.0;
                momentum += massI*molI->v();
                kE += 0.5*massI*magSqr(molI->v() - velocity);
                dof += molCloud_.cP().degreesOfFreedom(molI->id());

                const diagTensor& molMoI( molCloud_.cP().momentOfInertia(molI->id()));

                // angular speed
                const vector& molOmega(inv(molMoI) & molI->pi());
                angularKeSum += 0.5*(molOmega & molMoI & molOmega);

                stress_[cell] += massI*(molI->v() - velocity)*(molI->v() - velocity)
                                 + ((molOmega - angularVelocity) & molMoI) * (molOmega - angularVelocity)
                                    + 0.5*molI->rf(); // the virial is initially zero
            }
        }


        N_[cell] = mols;
        rhoN_[cell] = mols/volField[cell];
        rhoM_[cell] = mass/volField[cell];

        if(mass > 0.0)
        {
            USAM_[cell] = momentum/mass;
            UCAM_[cell] = momentum/mass;
        }

        if(dof > 0)
        {
            T_[cell] = (2.0*(kE+angularKeSum))/(kB*scalar(dof));
        }

        p_[cell] = (stress_[cell].xx()+stress_[cell].yy()+stress_[cell].zz())/(3.0*volField[cell]);
        stress_[cell] /= (volField[cell]);
    }

    N_.boundaryFieldRef() = N_.boundaryField().boundaryInternalField();
    rhoN_.boundaryFieldRef() = rhoN_.boundaryField().boundaryInternalField();
    rhoM_.boundaryFieldRef() = rhoM_.boundaryField().boundaryInternalField();
    USAM_.boundaryFieldRef() = USAM_.boundaryField().boundaryInternalField();
    UCAM_.boundaryFieldRef() = UCAM_.boundaryField().boundaryInternalField();
    T_.boundaryFieldRef() = T_.boundaryField().boundaryInternalField();
    p_.boundaryFieldRef() = p_.boundaryField().boundaryInternalField();
    stress_.boundaryFieldRef() = stress_.boundaryField().boundaryInternalField();

}


void polyPropertyFields::calculateField()
{
    molsV_ = 0.0;
    massV_ = 0.0;
    momV_ = vector::zero;
    velocity_ = vector::zero;
    angularSpeed_ = vector::zero;
    angularVelocity_ = vector::zero;

    forAll(molCloud_.cellOccupancy(), cell)
    {
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cell];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            if(findIndex(molIds_, molI->id()) != -1)
            {
//                 const polyMolecule::constantProperties& constProp = molCloud_.constProps(molI->id());
                const scalar& massI = molCloud_.cP().mass(molI->id());
                massV_[cell] += massI;
                momV_[cell] += massI*molI->v();

                const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                // angular speed
                const vector& molOmega(inv(molMoI) & molI->pi());

                angularSpeed_[cell] += molOmega;
            }
        }
    }


//     if(timeVel_.outputTime())
    {
        velocity_ = vector::zero;
        angularVelocity_ = vector::zero;

        forAll(velocity_, c)
        {
            if(massV_[c] > 0.0)
            {
               velocity_[c] = momV_[c]/massV_[c];
            }
            if(molsV_[c] > 0.0)
            {
                angularVelocity_[c] = angularSpeed_[c]/molsV_[c];
            }
        }
    }


    forAll(molCloud_.cellOccupancy(), cell)
    {
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cell];

        scalar mass = 0.0;
        vector momentum = vector::zero;

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            if(findIndex(molIds_, molI->id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(molI->id());

                mols_[cell] += 1.0;
                mass += massI;
                momentum += massI*molI->v();
                kE_[cell] += 0.5*massI*magSqr(molI->v() - velocity_[cell]);
                dof_[cell] += molCloud_.cP().degreesOfFreedom(molI->id());

                const diagTensor& molMoI( molCloud_.cP().momentOfInertia(molI->id()));
                const vector& molOmega(inv(molMoI) & molI->pi());
                angularKe_[cell] += 0.5*(molOmega & molMoI & molOmega);

                kineticTensor_[cell] += massI*(molI->v() - velocity_[cell])
                                                *(molI->v() - velocity_[cell])
                                        + (
                                            ( (molOmega - angularVelocity_[cell]) & molMoI)
                                            * (molOmega - angularVelocity_[cell])
                                            );

                virialTensor_[cell] += 0.5*molI->rf();
            }
        }

        mass_[cell] += mass;
        mom_[cell] += momentum;

        if(mass > 0.0)
        {
            velocityB_[cell] += momentum/mass;
        }
    }

    if(time_.outputTime())
    {
    	const scalar& nAvTimeSteps = nAvTimeSteps_;

        const scalar& kB = molCloud_.redUnits().kB();

        const scalarField& volField = mesh_.cellVolumes();

        forAll(mass_, cell)
        {
            N_[cell] = mols_[cell]/(nAvTimeSteps);
            rhoN_[cell] = mols_[cell]/(volField[cell]*nAvTimeSteps);
            rhoM_[cell] = mass_[cell]/(volField[cell]*nAvTimeSteps);
            USAM_[cell] = velocityB_[cell]/nAvTimeSteps;

            UCAM_[cell] = vector::zero;

            if(mass_[cell] > 0.0)
            {
                UCAM_[cell] = mom_[cell]/mass_[cell];
            }

            if(dof_[cell] > 0.0)
            {
                T_[cell] = (2.0*(kE_[cell]+angularKe_[cell]))/(kB*dof_[cell]);

                p_[cell] = tr( (3.0*mols_[cell]*kineticTensor_[cell]/dof_[cell]) + virialTensor_[cell])
                                        /( 3.0*volField[cell]*nAvTimeSteps );

                stress_[cell] = ((3.0*mols_[cell]*kineticTensor_[cell]/dof_[cell]) + virialTensor_[cell])
                                        /(volField[cell]*nAvTimeSteps);
            }
        }


        N_.boundaryFieldRef() = N_.boundaryField().boundaryInternalField();
        rhoN_.boundaryFieldRef() = rhoN_.boundaryField().boundaryInternalField();
        rhoM_.boundaryFieldRef() = rhoM_.boundaryField().boundaryInternalField();
        USAM_.boundaryFieldRef() = USAM_.boundaryField().boundaryInternalField();
        UCAM_.boundaryFieldRef() = UCAM_.boundaryField().boundaryInternalField();
        T_.boundaryFieldRef() = T_.boundaryField().boundaryInternalField();
        p_.boundaryFieldRef() = p_.boundaryField().boundaryInternalField();
        stress_.boundaryFieldRef() = stress_.boundaryField().boundaryInternalField();

        //- reset
        if(resetAtOutput_)
        {
            mols_ = 0.0;
            mass_ = 0.0;
            mom_ = vector::zero;
            velocityB_ = vector::zero;
            kE_ = 0.0;
            dof_ = 0.0;
            angularKe_ = 0.0;
            kineticTensor_ = tensor::zero;
            virialTensor_ = tensor::zero;
        }
    }
}


// const volScalarField& polyPropertyFields::densityField() const
// {
//     return rhoM_;
// }

void polyPropertyFields::writeField()
{}


/*
void polyPropertyFields::setVolumeMeasurements()
{
    if(coarseGrainedMeasurements_)
    {
        scalarField volumes(mesh_.nCells(), 0.0);

        //- sample measurement from real cells
        forAll(mesh_.cellVolumes(), cell)
        {
            volumes[cell] = mesh_.cellVolumes()[cell];
        }

        //- parallel-processing
        if(Pstream::parRun())
        {
            List<label> nCellsProcs(Pstream::nProcs(), 0);

            nCellsProcs[Pstream::myProcNo()] = mesh_.nCells();

            //- sending
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const label proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << mesh_.nCells();
                    }
                }
            }

            //- receiving
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    label nCellsProc;

                    const label proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> nCellsProc;
                    }

                    nCellsProcs[p] = nCellsProc;
                }
            }

            //- reconstruct the fields on all processors
            List<scalarField> volumesField(Pstream::nProcs());

            forAll(volumesField, p)
            {
                volumesField[p].setSize(nCellsProcs[p], 0.0);
            }

            forAll(volumes, c)
            {
                volumesField[Pstream::myProcNo()][c] = volumes[c];
            }

            //- sending
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const label proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << volumes;
                    }
                }
            }

            //- receiving
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField volumesProc;

                    const label proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> volumesProc;
                    }

                    forAll(volumesProc, c)
                    {
                        volumesField[p][c] = volumesProc[c];
                    }
                }
            }

            //- use referred cells to update from processor boundary cells
            const atomisticReferredCellList& referredCells = molCloud_.il().ril();

            forAll(referredCells, rIL)
            {
                const atomisticReferredCell& cellRef = referredCells[rIL];
                const label& procNo = cellRef.sourceProc();
                const labelList& cellIntList = cellRef.realCellsForInteraction();
                const label& sourceCell = cellRef.sourceCell();

                forAll(cellIntList, cIL)
                {
                    const label& realCellM = cellIntList[cIL];
                    volumes[realCellM] += volumesField[procNo][sourceCell];
                }
            }
        }
        else //- handle periodic boundaries
        {
            const atomisticReferredCellList& referredCells = molCloud_.il().ril();

            scalarField volumesPer(mesh_.nCells(), 0.0);

            forAll(referredCells, rIL)
            {
                const atomisticReferredCell& cellRef = referredCells[rIL];
                const labelList& cellIntList = cellRef.realCellsForInteraction();
                const label& sourceCell = cellRef.sourceCell();

                forAll(cellIntList, cIL)
                {
                    const label& realCellM = cellIntList[cIL];

                    volumesPer[realCellM] += volumes[sourceCell];
                }
            }

            volumes += volumesPer;
        }

        const List< DynamicList<atomisticMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

        forAll(cellOccupancy, cell)
        {
            const labelList& dICL = molCloud_.il().dil().fil()[cell];

            forAll(dICL, dCell)
            {
                const label cellJ = dICL[dCell];
                volumes[cell] += mesh_.cellVolumes()[cellJ];
            }
        }

        forAll(volumes_, c)
        {
            volumes_[c] = volumes[c];
        }
    }
    else
    {
        forAll(volumes_, c)
        {
            volumes_[c] = mesh_.cellVolumes()[c];
        }
    }
}


void polyPropertyFields::setCoarseGrainedMeasurements()
{
    const List< DynamicList<atomisticMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    scalarField mass(mesh_.nCells(), 0.0);

    //- sample measurement from real cells
    forAll(cellOccupancy, cell)
    {
        const List<atomisticMolecule*>& molsInCell = cellOccupancy[cell];

        forAll(molsInCell, mIC)
        {
            atomisticMolecule* molI = molsInCell[mIC];

            if((molI->id() == molId_) || !oneSpecie_)
            {
                const atomisticMolecule::constantProperties& constProp = molCloud_.constProps(molI->id());
                mass[cell] += constProp.mass();
            }
        }
    }

    //- parallel-processing
    if(Pstream::parRun())
    {
        List<label> nCellsProcs(Pstream::nProcs(), 0);

        nCellsProcs[Pstream::myProcNo()] = mesh_.nCells();

        //- sending
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const label proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << mesh_.nCells();
                }
            }
        }

        //- receiving
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label nCellsProc;

                const label proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> nCellsProc;
                }

                nCellsProcs[p] = nCellsProc;
            }
        }

        //- reconstruct the fields on all processors

        List<scalarField> massField(Pstream::nProcs());

        forAll(massField, p)
        {
            massField[p].setSize(nCellsProcs[p], 0.0);
        }

        forAll(mass, c)
        {
            massField[Pstream::myProcNo()][c] = mass[c];
        }

        //- sending
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const label proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << mass;
                }
            }
        }

        //- receiving
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalarField massProc;

                const label proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> massProc;
                }

                forAll(massProc, c)
                {
                    massField[p][c] = massProc[c];
                }
            }
        }


        //- use referred cells to update from processor boundary cells
        const atomisticReferredCellList& referredCells = molCloud_.il().ril();

        forAll(referredCells, rIL)
        {
            const atomisticReferredCell& cellRef = referredCells[rIL];
            const label& procNo = cellRef.sourceProc();
            const labelList& cellIntList = cellRef.realCellsForInteraction();
            const label& sourceCell = cellRef.sourceCell();

            forAll(cellIntList, cIL)
            {
                const label& realCellM = cellIntList[cIL];
                mass[realCellM] += massField[procNo][sourceCell];
            }
        }
    }
    else //- handle periodic boundaries
    {
        const atomisticReferredCellList& referredCells = molCloud_.il().ril();

        scalarField massPer(mesh_.nCells(), 0.0);

        forAll(referredCells, rIL)
        {
            const atomisticReferredCell& cellRef = referredCells[rIL];
            const labelList& cellIntList = cellRef.realCellsForInteraction();
            const label& sourceCell = cellRef.sourceCell();

            forAll(cellIntList, cIL)
            {
                const label& realCellM = cellIntList[cIL];

                massPer[realCellM] += mass[sourceCell];
            }
        }

        mass += massPer;
    }

    // - coarse graining involves also measurement from neighbouring cells,
    //   using the real interaction list (i.e. within ~rCut distance)

    forAll(cellOccupancy, cell)
    {
        const labelList& dICL = molCloud_.il().dil().fil()[cell];

        forAll(dICL, dCell)
        {
            const List< atomisticMolecule* >& molsInCellJ = cellOccupancy[dICL[dCell]];

            forAll(molsInCellJ, m)
            {
                atomisticMolecule* molJ = molsInCellJ[m];

                if((molJ->id() == molId_) || !oneSpecie_)
                {
                    const atomisticMolecule::constantProperties& constProp = molCloud_.constProps(molJ->id());
                    mass[cell] += constProp.mass();
                }
            }
        }
    }

    forAll(mass_, c)
    {
        mass_[c] += mass[c];
    }
}*/

void polyPropertyFields::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

// void polyPropertyFields::measureDuringForceComputation
// (
//     polyMolecule* molReal,
//     polyReferredMolecule* molRef
// ){}

void polyPropertyFields::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

// void polyPropertyFields::measureDuringForceComputationSite
// (
//     polyMolecule* molReal,
//     polyReferredMolecule* molRef,
//     label sReal,
//     label sRef
// ){}

const propertyField& polyPropertyFields::fields() const
{
    return fields_;
}

/*
void polyPropertyFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}
*/
} // End namespace Foam

// ************************************************************************* //
