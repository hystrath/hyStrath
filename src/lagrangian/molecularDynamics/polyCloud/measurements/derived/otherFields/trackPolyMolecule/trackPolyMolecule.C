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

#include "trackPolyMolecule.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(trackPolyMolecule, 0);

addToRunTimeSelectionTable(polyField, trackPolyMolecule, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
trackPolyMolecule::trackPolyMolecule
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
    molId_(-1),
    startPoint_(propsDict_.lookup("startPoint")),
    trackingNumber_(-1),
	trackingOrientation1_(2, vector::zero),
	trackingOrientation2_(2, vector::zero),
	trackingOrientation3_(2, vector::zero),
	trackingQ_(2, tensor::zero),
	trackingPi_(2, vector::zero),
	trackingLinMom_(2, vector::zero),
	trackingPositions_(2, vector::zero),
	trackingSitePositions_(2),
	phi_(2, 0.0),
	theta_(2, 0.0),
	cells_(2, -1),
    outputParaFoamFile_(false),
    outputVMDFile_(false)  
{
    const List<word>& idList(molCloud_.pot().idList());
    const word molId = propsDict_.lookup("molId");
    molId_ = findIndex(idList, molId);

    if(molId_ == -1)
    {
        FatalErrorIn("trackPolyMolecule::trackPolyMolecule()")
            << "Cannot find molId: " << molId << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    if (propsDict_.found("outputParaFoamFile"))
    {
        outputParaFoamFile_ = Switch(propsDict_.lookup("outputParaFoamFile"));
    }
    
    if (propsDict_.found("outputVMDFile"))
    {
        outputVMDFile_ = Switch(propsDict_.lookup("outputVMDFile"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

trackPolyMolecule::~trackPolyMolecule()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void trackPolyMolecule::createField()
{
    const polyMolecule::constantProperties& constProp = molCloud_.constProps(molId_);

    forAll(trackingSitePositions_, p)
    {
        trackingSitePositions_[p].setSize(constProp.sites().size(), vector::zero);
    }


    // find the closest polyMolecule to the given startPoint

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    label trackingNumber = -1;
    scalar shortestDistance = GREAT;
    vector startingMolPosition = vector::zero;

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        if(mol().id() == molId_)
        {
            const vector& rI = mol().position();
            scalar rSIMag = mag(rI - startPoint_);

            if(rSIMag < shortestDistance)
            {
                startingMolPosition = rI;
                shortestDistance = rSIMag;
                trackingNumber = mol().trackingNumber();
            }
        }
    }


    if(Pstream::parRun())
    {
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << trackingNumber << shortestDistance << startingMolPosition;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label trackingNumberProc;
                scalar shortestDistanceProc;
                vector startingMolPositionProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> trackingNumberProc >> shortestDistanceProc >> startingMolPositionProc;
                }

                if(shortestDistanceProc < shortestDistance)
                {
                    trackingNumber = trackingNumberProc;
                    shortestDistance = shortestDistanceProc;
                    startingMolPosition = startingMolPositionProc;
                }
                else if(shortestDistanceProc == shortestDistance)
                {
                    if(Pstream::myProcNo() > p)
                    {
                        trackingNumber = trackingNumberProc;
                        shortestDistance = shortestDistanceProc;
                        startingMolPosition = startingMolPositionProc;
                    }
                }
            }
        }
    }

//     cellI_ = mesh_.findCell(startingMolPosition);

    trackingNumber_ = trackingNumber;

    Info << "trackingMolecule: " << fieldName_ << ", starting position: " 
         << startingMolPosition << endl;
}


void trackPolyMolecule::calculateField()
{
    if(time_.outputTime())
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

	    tensor trackingQ = tensor::zero;
        vector trackingPi = vector::zero;
        vector trackingLinMom = vector::zero;		
        label cell = -1;
        vector trackingPosition = vector::zero;

        const polyMolecule::constantProperties& constProp = molCloud_.constProps(molId_);

        vectorField trackingMolSitePositions(constProp.sites().size(), vector::zero);

        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(mol().trackingNumber() == trackingNumber_)
            {
		        trackingQ = mol().Q();
		        trackingPi = mol().Q() & mol().pi();
		        const polyMolecule::constantProperties& constProp = molCloud_.constProps(mol().id());		
		        trackingLinMom = constProp.mass()*mol().v();		
                trackingPosition = mol().position();
                cell = mol().cell();

                forAll(trackingMolSitePositions, p)
                {
                    trackingMolSitePositions[p] = mol().sitePositions()[p];
                }
            }
        }


        forAll(trackingMolSitePositions, p)
        {
            trackingSitePositions_[1][p] = trackingMolSitePositions[p];
        }

        trackingQ_[1] = trackingQ;
        trackingPi_[1] = trackingPi;
        trackingLinMom_[1] = trackingLinMom;
        trackingPositions_[1] = trackingPosition;
        cells_[1] = cell;
       
        if(Pstream::parRun())
        {
            //- sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << trackingPosition << cell 
                                    << trackingMolSitePositions << trackingQ 
                                    << trackingPi << trackingLinMom;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
		            tensor trackingQProc;
		            vector trackingPiProc;
		            vector trackingLinMomProc;
                    vector trackingPositionProc;
                    label cellProc;
                    vectorField trackingMolSitePositionsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> trackingPositionProc >> cellProc 
                                      >> trackingMolSitePositionsProc >> trackingQProc 
                                      >> trackingPiProc >> trackingLinMomProc;
                    }

                    if(cellProc != -1)
                    {
			            trackingQ_[1] = trackingQProc;
			            trackingPi_[1] = trackingPiProc;
			            trackingLinMom_[1] = trackingLinMomProc;
                        trackingPositions_[1] = trackingPositionProc;
                        cells_[1] = cellProc;

                        forAll(trackingMolSitePositionsProc, p)
                        {
                            trackingSitePositions_[1][p] = trackingMolSitePositionsProc[p];
                        }
                    }
                }
            }
        }
    }
}



void trackPolyMolecule::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            fileName timePath(runTime.path()/runTime.timeName()/"uniform");

            writeTimeData
			(
				timePath,
				fieldName_+"_trackPositions",
				trackingPositions_,
				false
			);
			writeTimeData
			(
				timePath,
				fieldName_+"_trackingPi",
				trackingPi_,
				false
			);
			writeTimeData
			(
				timePath,
				fieldName_+"_trackingLinMom",
				trackingLinMom_,
				false
			);
	    
	        forAll(trackingQ_, p)
	        {
	            const tensor& q = trackingQ_[p];	
	        
	            trackingOrientation1_[p] = q & vector(1,0,0);
	            trackingOrientation2_[p] = q & vector(0,1,0);	    
	            trackingOrientation3_[p] = q & vector(0,0,1);    
	        }

	        writeTimeData
			(
				timePath,
				"trackOrientation2_"+fieldName_,
				trackingOrientation2_,
				false
			);
	    
	        //- a small script for visualising the trajectory (and other properties)
            //  of the polyMolecule in paraFOAM	
    
            if(outputParaFoamFile_)
            {
                OFstream positionsFile(timePath/fieldName_+"_positions.xyz");
                OFstream sitePositionsFile(timePath/fieldName_+"_sitePositions.xyz");
                OFstream siteIdsFile(timePath/fieldName_+"_siteIds.xy");
                OFstream bondsFile(timePath/fieldName_+"_bonds.xyz");
    
                OFstream QFile(timePath/fieldName_+"_Q.xyz");	
                OFstream piFile(timePath/fieldName_+"_pi.xyz");	
                OFstream linMomFile(timePath/fieldName_+"_linMom.xyz");	    
                OFstream timeFile(timePath/fieldName_+"_timeField.xy");
	            OFstream orientation1File(timePath/fieldName_+"_orientation1.xyz"); 
	            OFstream orientation2File(timePath/fieldName_+"_orientation2.xyz"); 
	            OFstream orientation3File(timePath/fieldName_+"_orientation3.xyz"); 	
                OFstream phiFile(timePath/fieldName_+"_phi.xy"); 
                OFstream thetaFile(timePath/fieldName_+"_theta.xy"); 
    
                label nPositions = 0;
    
                forAll(cells_, p)
                {
                    if(cells_[p] != -1)
                    {
                        nPositions++;
                    }
                }
    
                if (positionsFile.good() && QFile.good() && piFile.good() && linMomFile.good() 
	                && timeFile.good() /*&& orientation1File.good() */
	                && orientation2File.good() && phiFile.good() /*&& orientation3File.good()*/)
                {
                    positionsFile << nPositions << endl;
                    positionsFile << "(" << endl;
                    
                    const polyMolecule::constantProperties& constProp = molCloud_.constProps(molId_);
    
                    label nPositionsSites = label(nPositions*constProp.sites().size());
    
                    sitePositionsFile << nPositionsSites << endl;
                    sitePositionsFile << "(" << endl;
    
                    siteIdsFile << nPositionsSites << endl;
                    siteIdsFile << "(" << endl;
                    
                    bondsFile << nPositionsSites << endl;
                    bondsFile << "(" << endl;
    
                    QFile << nPositions << endl;
                    QFile << "(" << endl;		
    
                    piFile << nPositions << endl;
                    piFile << "(" << endl;			
		    
                    linMomFile << nPositions << endl;
                    linMomFile << "(" << endl;			
		    
                    orientation1File << nPositions << endl;
                    orientation1File << "(" << endl;		
		    
                    orientation2File << nPositions << endl;
                    orientation2File << "(" << endl;			
    
		            phiFile << nPositions << endl;
                    phiFile << "(" << endl;    
    
                    thetaFile << nPositions << endl;
                    thetaFile << "(" << endl;    
    
                    orientation3File << nPositions << endl;
                    orientation3File << "(" << endl;		
		    
                    timeFile << nPositions << endl;
                    timeFile << "(" << endl;
    
                    forAll(trackingPositions_, p)
                    {
                        const vector& r = trackingPositions_[p];
                        const vectorField& sitePos = trackingSitePositions_[p];
                        const tensor& q = trackingQ_[p];
                        const vector& lM = trackingLinMom_[p];
                        const vector& pi = trackingPi_[p];
                        const vector& tr1 = trackingOrientation1_[p];
                        const vector& tr = trackingOrientation2_[p];
                        const vector& tr3 = trackingOrientation3_[p];
    
                        if(cells_[p] != -1)
                        {
                            positionsFile 
                                << "(" << r.x() << " " << r.y() << " "
                                << r.z() << ") " << cells_[p]
                                << endl;
    
                            forAll(sitePos, s)
                            {
                                sitePositionsFile
                                    << "(" << sitePos[s].x() << " " << sitePos[s].y() << " "
                                    << sitePos[s].z() << ") " << cells_[p]
                                    << endl; 
    
                                siteIdsFile << constProp.sites()[s].siteId() << endl;
                            }
    
                            // bonds                                                         
                            bondsFile << "(" << (sitePos[2] - sitePos[0]).x() << " " 
                                    <<  (sitePos[2] - sitePos[0]).y() << " "
                                    << (sitePos[2] - sitePos[0]).z()
                                    << ") " << endl;
    
                            bondsFile << "(" << (sitePos[2] - sitePos[1]).x() << " " 
                                    <<  (sitePos[2] - sitePos[1]).y() << " "
                                    << (sitePos[2] - sitePos[1]).z()
                                    << ") " << endl;
    
                            bondsFile << "(0 0 0)" << endl;
                            bondsFile << "(0 0 0)" << endl;
    
                            QFile 
                                << "(" << q.xx() << " " << q.xy() << " " << q.xz() << " "
                                << q.yx() << " " << q.yy() << " " << q.yz() << " "
                                << q.zx() << " " << q.zy() << " " << q.zz() << " "
                                << ") " 
                                << endl;
    
                            piFile 
                                << "(" << pi.x() << " " << pi.y() << " "
                                << pi.z() << ") " 
                                << endl;
    
                            linMomFile 
                                << "(" << lM.x() << " " << lM.y() << " "
                                << lM.z() << ") " 
                                << endl;
    
                            orientation1File    
                                << "(" << tr1.x() << " " << tr1.y() << " "
                                << tr1.z() << ") " 
                                << endl;
                    
                            orientation2File    
                                << "(" << tr.x() << " " << tr.y() << " "
                                << tr.z() << ") " 
                                << endl;    
                    
                            orientation3File    
                                << "(" << tr3.x() << " " << tr3.y() << " "
                                << tr3.z() << ") "
                                << endl;          
    
                            scalar phi = acos(tr.z())*180.0/constant::mathematical::pi;
                            scalar theta = acos(tr.x()/sin(phi))*180.0/constant::mathematical::pi;
    
                            phiFile 
                                << phi << endl;
    
                            thetaFile 
                                << theta << endl;
                        }
                    }
    
                    positionsFile << ")" << endl;
                    sitePositionsFile << ")" << endl;
                    siteIdsFile << ")" << endl;
                    bondsFile << ")" << endl;
    
                    orientation1File << ")" << endl;    
                    orientation2File << ")" << endl;   		
                    orientation3File << ")" << endl;   	
		            QFile << ")" << endl;
		            piFile << ")" << endl;		
		            linMomFile << ")" << endl;			
		            phiFile << ")" << endl;
                    thetaFile << ")" << endl;
    
                }
                else
                {
                    FatalErrorIn( "trackPolyMolecule::writeField")
                        << "Cannot open file "
                        << positionsFile.name() << "or "<< QFile.name() << "or "<< piFile.name() << "or "<< linMomFile.name()
                        << "or "<< orientation2File.name()
                        << abort(FatalError);
                }
            }
            
	        if(outputVMDFile_)
            {
	        
	        }
        }
    }
}

void trackPolyMolecule::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void trackPolyMolecule::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& trackPolyMolecule::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
