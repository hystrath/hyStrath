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

#include "polyGrapheneSheetPeriodicSkew.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyGrapheneSheetPeriodicSkew, 0);

addToRunTimeSelectionTable(polyConfiguration, polyGrapheneSheetPeriodicSkew, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyGrapheneSheetPeriodicSkew::polyGrapheneSheetPeriodicSkew
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyConfiguration(molCloud, dict)
{
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyGrapheneSheetPeriodicSkew::~polyGrapheneSheetPeriodicSkew()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyGrapheneSheetPeriodicSkew::setInitialConfiguration()
{
    scalar b = readScalar(mdInitialiseDict_.lookup("bondLength"));
    
    const reducedUnits& rU = molCloud_.redUnits();
    
    b /= rU.refLength();
    
    scalar theta = 120*constant::mathematical::pi/180.0;
    scalar c = 2.0*b*sin(0.5*theta);

    Info << nl << "Information on properties of graphene sheet" << nl << endl;
    
    Info << "c = " << c
         << nl << "b = " << b
         << endl;    

    boundBox box = mesh_.bounds();
    
    vector nL = mdInitialiseDict_.lookup("lengthNormal");    
    vector nB = mdInitialiseDict_.lookup("breadthNormal");
    vector nZ = mdInitialiseDict_.lookup("heightNormal");    
    scalar Z = readScalar(mdInitialiseDict_.lookup("Z"));
    
    nL /= mag(nL);
    nB /= mag(nB);
    nZ /= mag(nZ);
    
    scalar Lx = box.span() & nL;
    scalar Ly = box.span() & nB;
    

    scalar Lscal = Lx/(3.0*b);
    scalar Bscal = Ly/c;
    
    Info<< "L (scalar) = " << Lscal
        << ", B (scalar) = " << Bscal 
        << endl;

    label L = label(Lscal+0.5);
    label B = label(Bscal+0.5);
    
    scalar b1 = Lx/(3*L);
    scalar b2 = Ly/(2*B*sin(0.5*theta) );
    scalar c2 = 2.0*b2*sin(0.5*theta);
    
    Info << "new b1 = " << b1
        << ", new b2 = " << b2
        << ", new c2 = " << c2
        << endl;

    
    // first layer of hexagon rings
    
    DynamicList<vector> Glayer;
    
    for (label i=1; i<L+1; i++)
    {
        vector r1 = (3*b1*i)*nL;
        vector r2 = (b1+(3*b1*i))*nL;
        vector r3 = (1.5*b1+(3*b1*i))*nL + 0.5*c2*nB;
        vector r4 = (b1+(3*b1*i))*nL + c2*nB;
        vector r5 = (3*b1*i)*nL + c2*nB;
        vector r6 = ((3*b1*i)-(0.5*b1))*nL + 0.5*c2*nB;
        
        Glayer.append(r1);
        Glayer.append(r2);
        Glayer.append(r3);
        Glayer.append(r4);
        Glayer.append(r5);
        Glayer.append(r6);
    }
    
    Glayer.shrink();
    
    // shift entire sheet so that the first atom which was created 
    // becomes centred at (0.0, 0.0)
    
    forAll (Glayer, a)
    {
       Glayer[a] -= nL*3*b1;
       Glayer[a] += nL*b1;
    }    
   
    // make a template for the first layer of hexagons
    DynamicList<vector> Gsheet;
    
    forAll (Glayer, a)
    {
        Gsheet.append(Glayer[a]);
    }
    
    // make a copy of the layer in the other direction to create the sheet
    
    for (label j=1; j<B; j++)
    {
        forAll (Glayer, a)
        {
            Gsheet.append(Glayer[a] + c2*j*nB);
        }
    }
    
    Gsheet.shrink();
    
    // remove any atoms that overlap eachother (to prevent the MD code from blowing up)
    
    DynamicList<vector> Gnew;
    
    scalar tolerance = 0.1;
    
    forAll (Gsheet, i)
    {
        const vector& rI = Gsheet[i];
        
        bool overlapping = false;
        
        forAll (Gnew, j)
        {
            const vector& rJ = Gnew[j];
            scalar rMag = mag(rI - rJ);
        
            if (rMag < tolerance)
            {
                overlapping = true;
            }
        }
    
        if (!overlapping)
        {
            Gnew.append(rI);
        }
    }
    
    Gnew.shrink();
    
    // shift sheet upwards so that there is a gap in the y direction
    
    forAll (Gnew, i)
    {
       Gnew[i] += nB*c2/4;
    }  

    // shift sheet in the z direction
    
    forAll (Gnew, i)
    {
        Gnew[i] += nZ*Z;
    }
    
    
    // READ IN MORE PROPERTIES FROM DICTIONARY

    word molIdName(mdInitialiseDict_.lookup("molId"));

    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polyGrapheneSheetPeriodicSkewId::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdName 
            << nl << "in moleculeProperties/idList."
            << exit(FatalError);
    }

    scalar temperature = 0.0;
    vector bulkVelocity = vector::zero;

    if (mdInitialiseDict_.found("temperature"))
    {
        temperature = readScalar(mdInitialiseDict_.lookup("temperature"));
    }

    if (mdInitialiseDict_.found("bulkVelocity"))
    {
        bulkVelocity = mdInitialiseDict_.lookup("bulkVelocity");
    }
    
    initialiseVelocities_ = false;
    
    if(mdInitialiseDict_.found("initialiseVelocities"))
    {
        initialiseVelocities_ = Switch(mdInitialiseDict_.lookup("initialiseVelocities"));
    }

    bool tethered = false;

    //- assume graphene sheet molecules are frozen always
    bool frozen = true;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }
    
    label noCatomsCreated = 0;
    
    forAll(Gnew, i)
    {  
        vector p = Gnew[i];
        
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            p,
            cell,
            tetFace,
            tetPt
        );

        if(cell != -1)
        {
            insertMolecule
            (
                p,
                cell,
                tetFace,
                tetPt,
                molId,
                tethered,
                frozen,
                temperature,
                bulkVelocity
            );
            noCatomsCreated++;
        }
        else
        {
            Info << "WARNING: molecule at position = " << p << ", out of mesh" << endl;
        }
    }

    Info << tab << "Graphene atoms created: " << noCatomsCreated << endl;
      
}


} // End namespace Foam

// ************************************************************************* //
