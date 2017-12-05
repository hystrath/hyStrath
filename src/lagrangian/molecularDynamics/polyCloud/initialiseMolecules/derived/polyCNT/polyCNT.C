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

#include "polyCNT.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyCNT, 0);

addToRunTimeSelectionTable(polyConfiguration, polyCNT, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCNT::polyCNT
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

polyCNT::~polyCNT()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyCNT::setInitialConfiguration()
{
    const reducedUnits& rU = molCloud_.redUnits();
    
    // READ IN PROPERTIES

    label N = readLabel(mdInitialiseDict_.lookup("N"));
    label M = readLabel(mdInitialiseDict_.lookup("M"));
    vector startPoint = mdInitialiseDict_.lookup("startPoint");
    vector endPoint = mdInitialiseDict_.lookup("endPoint");
    scalar bondLength = readScalar(mdInitialiseDict_.lookup("bondLengthSI"));
    
    bondLength /= rU.refLength();    
    
    vector cntNormal = (endPoint - startPoint) / mag(endPoint - startPoint);

    const scalar& bo = bondLength;

    scalar n = scalar(N);
    scalar m = scalar(M);

    word molIdName(mdInitialiseDict_.lookup("molId"));

    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polyCNTId::setInitialConfiguration()")
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

    bool tethered = false;

    bool stopBeforeCntRoll = false;

    //- ability to stop sheet before rolling, and view in its sheet format
    if (mdInitialiseDict_.found("stopCNTRoll"))
    {
        stopBeforeCntRoll = Switch(mdInitialiseDict_.lookup("stopCNTRoll"));
    }

    //- assume CNT molecules are frozen always
    bool frozen = true;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }

    // default value for triming the graphene sheet. Make sure it is a small number, 
    // smaller than the bond length

    scalar offset = 0.01;

    // over-right here to avoid compiling
    if (mdInitialiseDict_.found("trimSheetOffset"))
    {
        offset = readScalar(mdInitialiseDict_.lookup("trimSheetOffset"));
    }

    // you can specify a perpendicular vector 
    bool rotateAboutMidAxis = false;
    vector tPerp = vector::zero;

    if (mdInitialiseDict_.found("perpendicularVector"))
    {
        tPerp = mdInitialiseDict_.lookup("perpendicularVector");

        tPerp /= mag(tPerp);
        rotateAboutMidAxis = true;
        scalar rD = tPerp & cntNormal;

        if(rD > SMALL)
        {
            FatalErrorIn("void Foam::polyMoleculeCloud::createCNTs() : ")
                << "chosen perpendicularVector: " << tPerp 
                << " is not perpendicular to the cnt normal: " << cntNormal
                << ". Angle to normal is: " 
                << acos(rD)*180.0/constant::mathematical::pi 
                << " (degrees)."
                << nl
                << exit(FatalError);
        }
    }

    //Check that the chiral vector is correct
    if(N < M)
    {
        FatalErrorIn("void Foam::polyMoleculeCloud::createCNTs() : ")
            << " m is greater than n! "
            << exit(FatalError);
    }

    scalar length = mag(endPoint - startPoint);

    scalar theta = 
    atan
    (
        (sqrt(3.0)*n)/
        ((2.0*m)+ n)
    );

    vector a1 = vector(bo*sqrt(3.0), 0.0, 0.0);
    vector a2 = vector(0.5*bo*sqrt(3.0), 1.5*bo, 0.0);
    vector cH = n*a1 + m*a2;

//     scalar Ch = sqrt(3.0)*bondLength*sqrt(magSqr(n) + magSqr(m) + n*m);
//     scalar radius = (Ch*0.5)/constant::mathematical::pi;

    scalar D = mag(cH)/constant::mathematical::pi;

    // information 
    Info << nl << "Creating polyCNT. " << nl << endl;

    Info<< "CNT INFORMATION: "<< nl
        << "Radius: " << D*0.5 << nl
        << "Diameter: " << D << nl
        << "Length: " << length << nl
        << "Theta: " << theta << nl
        << "Chiral vector, cH: " << cH << nl
        << "a1: " << a1 << nl
        << "a2: " << a2 << nl
        << endl;


    // Create standard graphene sheet:

    scalar hexX = sqrt(3.0)*bo;

    scalar LxLeft = length*sin(theta) + hexX;
    scalar LxRight = mag(cH)*cos(theta) + hexX;
    scalar Ly = length*cos(theta) + mag(cH)*sin(theta) + 2.0*bo;

    DynamicList<vector> firstLayerPositions(0);

    label nColsRight = LxRight/hexX;
    label nColsLeft = LxLeft/hexX;
    label nRows = Ly/hexX;

    for(label i=0; i<nColsRight; i++)
    {
        vector rA = vector(hexX*i, 0.0, 0.0);
        vector rB = vector(hexX*(0.5+i), 0.5*bo, 0.0);
        vector rC = vector(hexX*(0.5+i), 1.5*bo, 0.0);
        vector rD = vector(hexX*i, 2.0*bo, 0.0);

        firstLayerPositions.append(rA);
        firstLayerPositions.append(rB);
        firstLayerPositions.append(rC);
        firstLayerPositions.append(rD);
    }

    for(label i=1; i<nColsLeft; i++)
    {
        vector rA = vector(hexX*i, 0.0, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);
        vector rB = vector(hexX*(0.5+i), 0.5*bo, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);
        vector rC = vector(hexX*(0.5+i), 1.5*bo, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);
        vector rD = vector(hexX*i, 2.0*bo, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);

        firstLayerPositions.append(rA);
        firstLayerPositions.append(rB);
        firstLayerPositions.append(rC);
        firstLayerPositions.append(rD);
    }

    DynamicList<vector> grapheneSheetPositions(0);

    forAll(firstLayerPositions, p)
    {
        grapheneSheetPositions.append(firstLayerPositions(p));
    }

    //firstLayerPositions.shrink();
    
    for(label i=1; i<nRows; i++)
    {
        forAll(firstLayerPositions, p)
        {
            grapheneSheetPositions.append(firstLayerPositions[p] + vector(0.0, 3.0*bo*i, 0.0));
        }
    }

    //grapheneSheetPositions.shrink();

    // rotate graphene sheet

    const tensor rotationMatrix
    (
        cos(theta), -sin(theta), 0.0,
        sin(theta), cos(theta), 0.0,
        0.0, 0.0, 1.0
    );

    vectorField rotatedGrapheneSheet;

    rotatedGrapheneSheet.transfer(grapheneSheetPositions);

    forAll(rotatedGrapheneSheet, p)
    {
        rotatedGrapheneSheet[p] = transform(rotationMatrix.T(), rotatedGrapheneSheet[p]);
    }

    // truncate graphene sheet

    vector pMin = vector::zero;
    vector pMax = vector(mag(cH), length, 0.0);

    // adding  buffer (to avoid truncation errors)
    vector span = pMax - pMin;
    span /= mag(span);

    pMin -= span*offset;
    pMax += span*offset;

    boundBox bb(pMin, pMax);

    DynamicList<vector> truncatedGrapheneSheet(0);

    forAll(rotatedGrapheneSheet, p)
    {
        if(bb.contains(rotatedGrapheneSheet[p]))
        {
            truncatedGrapheneSheet.append(rotatedGrapheneSheet[p]);
        }
    }

    //truncatedGrapheneSheet.shrink();

    if(stopBeforeCntRoll)
    {
        Info << "WARNING: CNT NOT ROLLED!" << endl;

        label noCatomsCreated = 0;
        
        forAll(truncatedGrapheneSheet, c)
        {  
            vector p = truncatedGrapheneSheet[c];
            
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

            noCatomsCreated += 1;
        }
    
        if (Pstream::parRun())
        {
            reduce(noCatomsCreated, sumOp<label>());
        }
    
        Info << tab << "CNT molecules added: " << noCatomsCreated << endl;
    }
    else
    {
        // Roll cnt
        vectorField rolledSheet;
        
        rolledSheet.transfer(truncatedGrapheneSheet);

        scalar s = mag(cH); // circumference
        scalar R = 0.5*D; // radius

        forAll(rolledSheet, c)
        {
            scalar phi = rolledSheet[c].x()*2.0*constant::mathematical::pi/s;

            vector pDash = vector(R*sin(phi), rolledSheet[c].y(), R*cos(phi));

            rolledSheet[c] = pDash;
        }

        // remove overlaps

        DynamicList<label> overlappingMolecules(0);

        forAll(rolledSheet, c1)
        {
            const vector& p1 = rolledSheet[c1];

            forAll(rolledSheet, c2)
            {
                if(c2 > c1)
                {
                    const vector& p2 = rolledSheet[c2];
                    
                    scalar rD = mag(p2-p1);

                    if(rD < 0.5*bo)
                    {
                        if(findIndex(overlappingMolecules, c2) == -1)
                        {
                            overlappingMolecules.append(c2);
                        }
                    }
                }
            }
        }

        //overlappingMolecules.shrink();

        label Ncnt = (rolledSheet.size() - overlappingMolecules.size());

        vectorField cntMolecules(Ncnt, vector::zero);

        label counter = 0;

        forAll(rolledSheet, c)
        {
            if(findIndex(overlappingMolecules, c) == -1)
            {
                cntMolecules[counter] = rolledSheet[c];
                counter++;
            }
        }

        // orient CNT from the local co-ordinate axis used to setup CNT, to fit 
        // in the global co-ordinate axis where it was defined.


        vector t1 = vector::zero;
        vector t2 = vector::zero;

        if(rotateAboutMidAxis)
        {
            t1 = tPerp;
            t2 = cntNormal ^ t1;
        }
        else
        {
            scalar magV = 0.0;
            vector tangent;
        
            while (magV < SMALL)
            {
                vector testThis = molCloud_.rndGen().sampleVectorMD<vector>();
        
                tangent = testThis - (testThis & cntNormal)*cntNormal;
                magV = mag(tangent);
            }
    
            t1 = tangent/magV;
            t2 = cntNormal ^ t1;
        }

        Info << "Tangential vectors - t1: " << t1 << ", t2 : " << t2 << endl;

        forAll(cntMolecules, c)
        {
            cntMolecules[c] = startPoint + cntMolecules[c].y()*cntNormal 
                                + t1*cntMolecules[c].x() 
                                + t2*cntMolecules[c].z();
        }

        // create atoms

        label noCatomsCreated = 0;
        
        forAll(cntMolecules, c)
        {
            vector p = cntMolecules[c];

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
        
            if (cell != -1)
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

                noCatomsCreated += 1;
            }
            else
            {
                Info << "WARNING - Atom at specified position " << p
                     << ", does not correspond to a mesh cell.... deleting"
                     << nl
                     << endl;
            }
        }
    
        if (Pstream::parRun())
        {
            reduce(noCatomsCreated, sumOp<label>());
        }
    
        Info << tab << " CNT molecules added: " << noCatomsCreated << endl;
    }
}




} // End namespace Foam

// ************************************************************************* //
