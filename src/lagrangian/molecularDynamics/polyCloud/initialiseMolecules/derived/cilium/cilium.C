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

#include "cilium.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(cilium, 0);

addToRunTimeSelectionTable(polyConfiguration, cilium, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cilium::cilium
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyConfiguration(molCloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cilium::~cilium()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void cilium::setInitialConfiguration()
{
    label initialSize = molCloud_.size();

    Info << nl << "Creating cilium " << endl;

    const word molIdName(mdInitialiseDict_.lookup("molId")); 
    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("cilium::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdName << nl << "in idList."
            << exit(FatalError);
    }

    const reducedUnits& rU = molCloud_.redUnits();
    
    scalar temperature = 300/rU.refTemp();
    vector bulkVelocity = vector::zero;
    
    //- start point is the fixed point
    vector startPoint = mdInitialiseDict_.lookup("startPoint");
    
    vector endPoint = mdInitialiseDict_.lookup("endPoint");
    

    
    
    bool tethered = false;
    
    scalar spacing = readScalar(mdInitialiseDict_.lookup("spacingSI"));
    
    spacing /= rU.refLength();

    DynamicList<vector> positions;
    
    // unit vector 
    vector n = endPoint - startPoint;
    
    scalar magSE = mag(n);
    
    n /= mag(n);
    
    scalar Ns = magSE/spacing;
    
    label N = label(Ns) + 1;
    
    // new spacing
    scalar s = magSE/scalar(N-1);
    
    for (label i = 0; i < N; i++)
    {
        vector p = startPoint + i*s*n;
        positions.append(p);
    }
    
    positions.shrink();

    Info << nl << " No of sites found = " << positions.size() << endl;

    DynamicList<label> frozenAtoms(0);
    
    if(mdInitialiseDict_.found("frozenAtoms"))
    {
        List<label> molecules = List<label>(mdInitialiseDict_.lookup("frozenAtoms"));
        
        if(molecules.size() > positions.size())
        {
            FatalErrorIn("cilium::setInitialConfiguration()")
                << "You can't have more frozen atoms than you have atoms = " << positions.size() 
                << exit(FatalError);            
        }
        
        forAll(molecules, i)
        {
            frozenAtoms.append(molecules[i]);
        }
        
        Info << "frozen atoms = " << frozenAtoms << endl;
    }
    
    
    
    
    // insert molecules in cloud
    label nMolsInserted = 0;
    
    forAll(positions, i)
    {
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;
        
        mesh_.findCellFacePt
        (
            positions[i],
            cell,
            tetFace,
            tetPt
        );
        
        bool frozen = false;
        
        if(findIndex(frozenAtoms, i) != -1)
        {
            frozen = true;
        }
        
        if(cell != -1)
        {
            insertMoleculeLocal
            (
                positions[i],
                cell,
                tetFace,
                tetPt,
                molId,
                tethered,
                frozen,
                temperature,
                bulkVelocity
            );
            
            nMolsInserted++;
        }
    }

    Info<< nl << " No of initial cloud = " << initialSize
        << ", no of molecules inserted = " << nMolsInserted
        << endl;
    
        
    // write out of ordered locations 
    
    forAll(trackingNumbers_, i)
    {
        if(i < trackingNumbers_.size() - 1)
        {
            trackingNumbersA_.append(trackingNumbers_[i]);
            trackingNumbersB_.append(trackingNumbers_[i+1]);
        }
    }

    Info << "trackingNumbersA = " << trackingNumbersA_ << endl;
    Info << "trackingNumbersB = " << trackingNumbersB_ << endl;
}


// simple modifications to the insert molecule function
void cilium::insertMoleculeLocal
(
    const point& position,
    const label cell,
    const label tetFace,
    const label tetPt, 
    const label& id, 
    const bool& tethered,
    const bool& frozen,
    const scalar& temperature,
    const vector& bulkVelocity
)
{
    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_TETHERED;
    }

    if (frozen)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_FROZEN;
    }

//     const polyMolecule::constantProperties& cP = molCloud_.constProps(id);

//     vector v = equipartitionLinearVelocity(temperature, molCloud_.cP().mass(id));

//     v += bulkVelocity;
    
    vector v = bulkVelocity;
    
    vector pi = vector::zero;

    tensor Q = I;

    if (!molCloud_.cP().pointMolecule(id))
    {
//         Info << "temperature = " << temperature << ", id = " << id << endl;
        
        pi = equipartitionAngularMomentum(temperature, id);
        scalar phi(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
        scalar theta(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
        scalar psi(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);

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
    }
    
    label tNI = molCloud_.getTrackingNumber();
    trackingNumbers_.append(tNI);
//     positions_.append(position);
    
    molCloud_.createMolecule
    (
        position,
        cell,
        tetFace,
        tetPt,     
        Q,
        v,
        vector::zero,
        pi,
        vector::zero,
        specialPosition,
        special,
        id,
        1.0,
        tNI
    );
    
}




} // End namespace Foam

// ************************************************************************* //
