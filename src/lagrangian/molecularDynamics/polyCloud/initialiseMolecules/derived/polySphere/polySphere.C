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

#include "polySphere.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polySphere, 0);

addToRunTimeSelectionTable(polyConfiguration, polySphere, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polySphere::polySphere
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polySphere::~polySphere()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polySphere::setInitialConfiguration()
{
    label initialSize = molCloud_.size();

    Info << nl << "Building quick simple lattice: " << endl;

    const scalar temperature(readScalar(mdInitialiseDict_.lookup("temperature")));

    const vector bulkVelocity(mdInitialiseDict_.lookup("bulkVelocity"));
    
    const vector midPoint(mdInitialiseDict_.lookup("midPoint"));
    
    const scalar radius(readScalar(mdInitialiseDict_.lookup("radius")));

    bool frozen = false;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }

    bool tethered = false;

    if (mdInitialiseDict_.found("tethered"))
    {
        tethered = Switch(mdInitialiseDict_.lookup("tethered"));
    }    
    
    const word molIdName(mdInitialiseDict_.lookup("molId")); 
    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polySphere::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdName << nl << "in idList."
            << exit(FatalError);
    }
    
    bool hemisphere = false;
    vector hemisphereVector = vector::zero;
    
    if (mdInitialiseDict_.found("hemisphereVector"))    
    {
        hemisphere = true;
        hemisphereVector = mdInitialiseDict_.lookup("hemisphereVector");        
        hemisphereVector /= mag(hemisphereVector);        
    }
    

    
    
    scalar numberDensity = 0.0;
        
    if (mdInitialiseDict_.found("massDensity"))
    {
//         const polyMolecule::constantProperties& cP(molCloud_.constProps(molId));

        scalar mass = molCloud_.cP().mass(molId);

        Info << "mass: " << mass << endl;

        scalar massDensity = readScalar
        (
            mdInitialiseDict_.lookup("massDensity")
        );

        numberDensity = massDensity / mass;

        if (massDensity < VSMALL)
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "massDensity too small, not filling zone "
                << abort(FatalError);
        }
    }
    else if (mdInitialiseDict_.found("massDensitySI"))
    {
//         const polyMolecule::constantProperties& cP(molCloud_.constProps(molId));

        scalar mass = molCloud_.cP().mass(molId);

        Info << "mass: " << mass << endl;

        scalar massDensity = readScalar
        (
            mdInitialiseDict_.lookup("massDensitySI")
        );
        
        const reducedUnits& rU =molCloud_.redUnits();

        massDensity /= rU.refMassDensity();
        
        numberDensity = massDensity / mass;
        
        Info << " number density in reduced units = " << numberDensity << endl;
        
        if (massDensity < VSMALL)
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "massDensity too small, not filling zone "
                << abort(FatalError);
        }
    }    
    

    
    scalar spacing = pow(  (1.0/numberDensity), (1.0/3.0) );
    
    vector xVector = vector(1, 0, 0);
    vector yVector = vector(0, 1, 0);
    vector zVector = vector(0, 0, 1);    

    label nMols = label(2*radius/spacing) + 1;
    
    vector startPoint = midPoint - xVector*radius - yVector*radius - zVector*radius;
        
    vector globalPosition = vector::zero;
    
    
    scalar x;
    scalar y;
    scalar z;
    
    for (label i = 0; i < nMols; i++)
    {
        x = i*spacing;

        for (label j = 0; j < nMols; j++)
        {
            y = j*spacing;

            for (label k = 0; k < nMols; k++)
            {
                z = k*spacing;

                globalPosition = startPoint 
                                + xVector*x 
                                + yVector*y 
                                + zVector*z;
                                
                scalar rIM = mag(globalPosition - midPoint);
                
                if(rIM <= radius)
                {
                    if(!hemisphere)
                    {
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
                                molId,
                                tethered,
                                frozen,
                                temperature,
                                bulkVelocity
                            );
                            
    //                         nMolsAdded++
                        }
                    }
                    else if(((globalPosition - midPoint) & hemisphereVector) > 0)
                    {
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
                                molId,
                                tethered,
                                frozen,
                                temperature,
                                bulkVelocity
                            );
                            
    //                         nMolsAdded++
                        }                    
                    }
                }
            }
        }
    }

    label finalSize = molCloud_.size();

    label nMolsAdded = finalSize - initialSize;

    if (Pstream::parRun())
    {
        reduce(nMolsAdded, sumOp<label>());
    }

    Info << tab << " molecules added: " << nMolsAdded << endl;

 
}



} // End namespace Foam

// ************************************************************************* //
