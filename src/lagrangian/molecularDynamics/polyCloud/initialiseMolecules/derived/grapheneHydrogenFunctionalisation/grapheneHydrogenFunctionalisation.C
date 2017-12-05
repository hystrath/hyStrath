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

#include "grapheneHydrogenFunctionalisation.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(grapheneHydrogenFunctionalisation, 0);

addToRunTimeSelectionTable(polyConfiguration, grapheneHydrogenFunctionalisation, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
grapheneHydrogenFunctionalisation::grapheneHydrogenFunctionalisation
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{

    
    

    option_ = "fixedPropertiesFromDict";
    
    if (mdInitialiseDict_.found("option"))
    {    
        const word option(mdInitialiseDict_.lookup("option"));    
        option_ = option;
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

grapheneHydrogenFunctionalisation::~grapheneHydrogenFunctionalisation()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void grapheneHydrogenFunctionalisation::setInitialConfiguration()
{
    label initialSize = molCloud_.size();
    
    if(option_ == "fixedPropertiesFromDict")
    {    
        fixedPropertiesFromDict();
    }
    
    label finalSize = molCloud_.size();

    nMolsAdded_ = finalSize - initialSize;

    if (Pstream::parRun())
    {
        reduce(nMolsAdded_, sumOp<label>());
    }

    Info << tab << " molecules added: " << nMolsAdded_ << endl; 
}

void grapheneHydrogenFunctionalisation::fixedPropertiesFromDict()
{

    // points on graphene 
    List<vector> molPoints = List<vector>(mdInitialiseDict_.lookup("molPoints"));
    
    // check for overlaps
    scalar threshold = 0.01;
    
    forAll(molPoints, i)
    {
        forAll(molPoints, j)
        {
            if(i != j)
            {
                scalar magRIJ = mag(molPoints[i]-molPoints[j]);
                
                if(magRIJ < threshold)
                {
                    FatalErrorIn("grapheneHydrogenFunctionalisation::setInitialConfiguration()")
                        << "You've specified similar atoms in the mdInitialiseDict: " 
                        << molPoints[i]
                        << ", " << molPoints[j]
                        << ", magRIJ = " << magRIJ
                        << exit(FatalError);
                }
            }
        }
    }
    
    
//     Info << "molPoints" << molPoints << endl;
    
    label nMols = molPoints.size();

    List<label> molIds(nMols, 0);
    List<bool> tetheredMols(nMols, false);
    List<bool> frozenMols(nMols, false);
//     List<scalar> temperatureMols(nMols, 0.0);
    List<vector> velocityMols(nMols, vector::zero);
    List<scalar> phiMols(nMols, 0.0);
    List<scalar> thetaMols(nMols, 0.0);
    List<scalar> psiMols(nMols, 0.0);


    {
        word molIdName(mdInitialiseDict_.lookup("molIdCarbon"));

        const List<word>& idList(molCloud_.cP().molIds());

        label molId = findIndex(idList, molIdName);

        if(molId == -1)
        {
            FatalErrorIn("grapheneHydrogenFunctionalisation::setInitialConfiguration()")
                << "Cannot find molecule id: " << molIdName 
                << nl << "in moleculeProperties/idList."
                << exit(FatalError);
        }
        
        molIdC_ = molId;
    }
    {
        word molIdName(mdInitialiseDict_.lookup("molIdCarbonNew"));

        const List<word>& idList(molCloud_.cP().molIds());

        label molId = findIndex(idList, molIdName);

        if(molId == -1)
        {
            FatalErrorIn("grapheneHydrogenFunctionalisation::setInitialConfiguration()")
                << "Cannot find molecule id: " << molIdName 
                << nl << "in moleculeProperties/idList."
                << exit(FatalError);
        }
        
        molIdCN_ = molId;
    }    
    {
        word molIdName(mdInitialiseDict_.lookup("molIdHydrogen"));

        const List<word>& idList(molCloud_.cP().molIds());

        label molId = findIndex(idList, molIdName);

        if(molId == -1)
        {
            FatalErrorIn("grapheneHydrogenFunctionalisation::setInitialConfiguration()")
                << "Cannot find molecule id: " << molIdName 
                << nl << "in moleculeProperties/idList."
                << exit(FatalError);
        }
        
        molIdH_ = molId;
    }    

    dr_ = readScalar(mdInitialiseDict_.lookup("dr_SI"));    
    
    const reducedUnits& rU = molCloud_.redUnits();
    
    dr_ /= rU.refLength();
    
//         const scalar T(readScalar(mdInitialiseDict_.lookup("temperature")));
//     const vector U(mdInitialiseDict_.lookup("velocity"));
    
    vector U = vector::zero;

    if (mdInitialiseDict_.found("velocity"))
    {
        vector vel(mdInitialiseDict_.lookup("velocity"));
        U = vel;
    }    
    
    bool frozen = true;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }

    bool tethered = false;

    if (mdInitialiseDict_.found("tethered"))
    {
        tethered = Switch(mdInitialiseDict_.lookup("tethered"));
    }

    scalar phi = 0.0;

    if (mdInitialiseDict_.found("phi"))
    {
        phi = readScalar(mdInitialiseDict_.lookup("phi"));
    }

    scalar theta = 0.0;

    if (mdInitialiseDict_.found("theta"))
    {
        theta = readScalar(mdInitialiseDict_.lookup("theta"));
    }

    scalar psi = 0.0;

    if (mdInitialiseDict_.found("psi"))
    {
        psi = readScalar(mdInitialiseDict_.lookup("psi"));
    }

    
    label nHydrogenInserted = 0;    
    
    // For each carbon find it, and also find it's two neighbours 
    // and then insert a hydrogen ion dr away 
    
    DynamicList<polyMolecule*> changeC;
    
    
    forAll(molPoints, i)
    {    
        DynamicList<polyMolecule*> molsA;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
       
        scalar rMin = GREAT;
        
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(mol().id() == molIdC_)
            {
                vector rT = molPoints[i];
                scalar magRIJ = mag(rT-mol().position());
                
                if(magRIJ < rMin)
                {
                    molsA.clear();
                    rMin = magRIJ;
                    polyMolecule* molI = &mol();
                    molsA.append(molI);
                }
            }
        }
        
        rMin = GREAT;        
        DynamicList<polyMolecule*> molsB;
        
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(mol().id() == molIdC_)
            {
                vector rT = molPoints[i];
                scalar magRIJ = mag(rT-mol().position());
                
                if(magRIJ < rMin)
                {
                    if(mol().trackingNumber() != molsA[0]->trackingNumber())
                    {
                        if(molsB.size() > 0)
                        {

                                molsB.clear();
                                rMin = magRIJ;
                                polyMolecule* molI = &mol();
                                molsB.append(molI);
                        }
                        else
                        {
                            molsB.clear();
                            rMin = magRIJ;
                            polyMolecule* molI = &mol();
                            molsB.append(molI);
                        }
                    }
                }
            }
        }        

        rMin = GREAT;        
        DynamicList<polyMolecule*> molsC;
        
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(mol().id() == molIdC_)
            {
                vector rT = molPoints[i];
                scalar magRIJ = mag(rT-mol().position());
                
                if(magRIJ < rMin)
                {
                    if
                    (
                        (mol().trackingNumber() != molsA[0]->trackingNumber()) &&
                        (mol().trackingNumber() != molsB[0]->trackingNumber())
                    )
                    {
                        if(molsC.size() > 0)
                        {
                            molsC.clear();
                            rMin = magRIJ;
                            polyMolecule* molI = &mol();
                            molsC.append(molI);
                        }
                        else
                        {
                            molsC.clear();
                            rMin = magRIJ;
                            polyMolecule* molI = &mol();
                            molsC.append(molI);
                        }
                    }
                }
            }
        }        

        vector rA = molsA[0]->position();
        vector rB = molsB[0]->position();
        vector rC = molsC[0]->position();
        
//         Info << "target position = " << molPoints[i]
//             << ", molsA = " << molsA[0]->position()
//             << ", molsB = " << molsB[0]->position()
//             << ", molsC = " << molsC[0]->position()
//             << " rAB = " << mag(molsA[0]->position() -  molsB[0]->position())
//             << " rAC = " << mag(molsA[0]->position() -  molsC[0]->position())
//             << endl;

        // average of all positions
        vector rM = (rA + rB + rC) / 3;
        
        vector unit = rA - rM;
        unit /= mag(unit);
        
        // switch mol A to molIdCNew_
//         molsA[0]->id() = molIdCN_;
        changeC.append(molsA[0]);
                
        // insert hydrogen atom
        
        vector globalPosition = rA + unit*dr_;
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            globalPosition,
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            insertMolecule
            (
                globalPosition,
                cell,
                tetFace,
                tetPt,
                molIdH_,
                tethered,
                frozen,
                phi*constant::mathematical::pi/180.0,
                theta*constant::mathematical::pi/180.0,
                psi*constant::mathematical::pi/180.0,
//                 temperatureMols[i],
                U
            );
            
            nHydrogenInserted++;
        }
        else
        {
            FatalErrorIn("Foam::grapheneHydrogenFunctionalisation::setInitialConfiguration()")
                << "Molecule position: " << globalPosition 
                << " is not located in the mesh." << nl
                << abort(FatalError);
        }        
    }
    
    forAll(changeC, i)
    {
        changeC[i]->id() = molIdCN_;
    }
    
    
}


void grapheneHydrogenFunctionalisation::insertMolecule
(
    const point& position,
    const label& cell,
    const label& tetFace,
    const label& tetPt,
    const label& id,
    const bool& tethered,
    const bool& frozen,
    const scalar& phi,
    const scalar& theta,
    const scalar& psi,
//     const scalar& temperature,
    const vector& velocity
)
{
    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_TETHERED;
    }

    if (frozen)//****
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_FROZEN;
    }

//     const polyMolecule::constantProperties& cP = molCloud_.constProps(id);

//     vector v = equipartitionLinearVelocity(temperature, cP.mass());
// 
//     v += bulkVelocity;

    vector pi = vector::zero;

    tensor Q = I;

//     scalar phi = 0.5*constant::mathematical::pi;
//     scalar theta = 0.0;
//     scalar psi = 0.0;

/*
    scalar phi(rndGen_.scalar01()*constant::mathematical::twoPi);

    scalar theta(rndGen_.scalar01()*constant::mathematical::twoPi);

    scalar psi(rndGen_.scalar01()*constant::mathematical::twoPi);*/

//     Info << "phi:" << phi << endl;

    Q = tensor
    (
        cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
        cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
        sin(psi)*sin(theta),
        - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
        - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
        cos(psi)*sin(theta),
        sin(theta)*sin(phi),
        - sin(theta)*cos(phi),
        cos(theta)
    );

    molCloud_.createMolecule
    (
        position,
        cell,
        tetFace,
        tetPt,
        Q,
        velocity,
        vector::zero,
        pi,
        vector::zero,
        specialPosition,
        special,
        id,
        1.0,
        molCloud_.getTrackingNumber()
    );
}

} // End namespace Foam

// ************************************************************************* //
