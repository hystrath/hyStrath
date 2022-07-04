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

Class
    Foam::cellInteractions<ParticleType>

Description

\*----------------------------------------------------------------------------*/

#include "cellInteractions.H"
#include "polyBoundaryMeshEntries.H"
#include "processorCyclicPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh

template<class ParticleType>
Foam::cellInteractions<ParticleType>::cellInteractions
(
    const polyMesh& mesh,
    const reducedUnits& rU,
    const cyclicBoundaries& cyclics,
    const scalar& rCut,
    const word& fieldName
)
:
    mesh_(mesh),
    cyclics_(cyclics),
    redUnits_(rU),
    fieldName_(fieldName),
    cloud_(mesh_, "referredParticleCloud", IDLList<ParticleType>()),
    optimised_(true),
    rCut_(rCut),
    dil_(mesh_.nCells()),
    fil_(mesh_.nCells()),
    inverseDIL_(mesh_.nCells()),
    recRefIds_(Pstream::nProcs()),
    sendSrcCells_(Pstream::nProcs()),
    refCells_(),
    sourceCellToRefs_(mesh_.nCells()),
    inverseFRIL_(mesh_.nCells()),
    ripl_(mesh_.nCells()),
//     ipl_(mesh_.nCells()),
    referredCloud_(),
    nCoupledPatches_(0),
    nCyclicProcessorBoundaries_(0),
    nProcessorBoundaries_(0),
    nCyclicBoundaries_(0),
    write_(false),
    writeReferredCells_(false)
{
    IOdictionary propDict
    (
        IOobject
        (
            "potentialDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (propDict.found("optimisedMesh"))
    {
        optimised_ = Switch(propDict.lookup("optimisedMesh"));
    }

    if (propDict.found("writeReferredCells"))
    {
        writeReferredCells_ = Switch(propDict.lookup("writeReferredCells"));
    }

    if (propDict.found("writeReferredCloud"))
    {
        write_ = Switch(propDict.lookup("writeReferredCloud"));
    }

    if(optimised_)
    {
        Info << nl << "Mesh: optimised" << nl << endl;
        checkMesh(rCut);
    }
    else
    {
        Info << nl << "Mesh: non-optimised" << nl << endl;
    }

    buildInteractionLists();
}


template<class ParticleType>
Foam::cellInteractions<ParticleType>::cellInteractions
(
    const polyMesh& mesh,
    const reducedUnits& rU,
    const cyclicBoundaries& cyclics,
    const label& dummy
)
:
    mesh_(mesh),
    redUnits_(rU),
    cyclics_(cyclics),
    cloud_(mesh_, "referredParticleCloud", IDLList<ParticleType>()),
    write_(false),
    writeReferredCells_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::cellInteractions<ParticleType>::~cellInteractions()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::cellInteractions<ParticleType>::buildInteractionLists()
{
    // build direct interaction list

    if(optimised_)
    {
        buildOptimisedDIL();
    }
    else
    {
        buildComplexDIL();
    }

    buildInverseDirectInteractionList();

    buildFullInteractionList();

    buildReferredCells();
}


template<class ParticleType>
void Foam::cellInteractions<ParticleType>::buildOptimisedDIL()
{
    Info << "Building direct interaction lists  - Your mesh is optimised for MD ..." << endl;

    List<DynamicList<label> > dil(mesh_.nCells());

    forAll(dil, c)
    {
        labelList cells = returnNeighbouringCells(c);

        forAll(cells, cN)
        {
            label cellN = cells[cN];

            if(includeCell(dil, c, cellN))
            {
                dil[c].append(cellN);
            }
        }
    }

    forAll(dil, c)
    {
        dil[c].shrink();
        dil_[c].setSize(dil[c].size());
        dil_[c].transfer(dil[c]);
    }

    // TESTING

//     Info << "dil = " << dil_ << endl;
//
//     forAll(dil_, c)
//     {
//         Info << "Cell Centre = " << mesh_.cellCentres()[c] << endl;
//         forAll(dil_[c], nC)
//         {
//             const label& cellN = dil_[c][nC];
//             Info << "   n = " << mesh_.cellCentres()[cellN]
//                 << endl;
//         }
//     }

    // TESTING


    Info << " ... done." << endl;
}


template<class ParticleType>
void Foam::cellInteractions<ParticleType>::buildComplexDIL()
{
    // to make sure the user sees this
    for (label j = 0; j < 50; j++)
    {
        Info << "WARNING - YOUR MESH IS NOT OPTIMISED FOR MD..." << endl;
    }

    Info << nl << "Building direct interacion lists (arbitrary mesh version)" << endl;

    List<DynamicList<label> > dil(mesh_.nCells());

    forAll(dil, c1)
    {
        forAll(dil, c2)
        {
            if(c1 > c2)
            {
                if(interactingCells(c1, c2, rCut_) )
                {
                    if(includeCell(dil, c1, c2) )
                    {
                        dil[c1].append(c2);
                    }
                }
            }
        }
    }

    forAll(dil, c)
    {
        dil[c].shrink();
        dil_[c].setSize(dil[c].size());
        dil_[c].transfer(dil[c]);
    }

    Info << " ... done." << endl;
}

// build full interaction list
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::buildFullInteractionList()
{
    List<DynamicList<label> > fil(mesh_.nCells());

    forAll(dil_, cellA)
    {
        const labelList& cellList = dil_[cellA];

        forAll(cellList, cell)
        {
            const label& cellB = cellList[cell];

            fil[cellA].append(cellB);
            fil[cellB].append(cellA);
        }
    }

    forAll(fil, c)
    {
        fil[c].shrink();
        fil_[c].setSize(fil[c].size());
        fil_[c].transfer(fil[c]);
    }
}

// build inverse direct interaction list
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::buildInverseDirectInteractionList()
{
    List<DynamicList<label> > inverseDIL(mesh_.nCells());

    forAll(dil_, d)
    {
        forAll(dil_[d], iC)
        {
            const label& cellJ = dil_[d][iC];
            inverseDIL[cellJ].append(d);
        }
    }

    forAll(inverseDIL, d)
    {
        inverseDIL[d].shrink();
        inverseDIL_[d].setSize(inverseDIL[d].size());
        inverseDIL_[d].transfer(inverseDIL[d]);
    }

//         Info << " inverse DIL: " << inverseDIL_ << endl;
}


template<class ParticleType>
Foam::labelList Foam::cellInteractions<ParticleType>::returnNeighbouringCells(const label& cellI)
{
//     Info << "cellI: " << cellI
//         << ", cell-centre: "
//         << mesh_.cellCentres()[cellI]
//         << endl;

    DynamicList<label> neighbCells(0);

    const labelList& points = mesh_.cellPoints()[cellI];

    forAll(points, p)
    {
        const label& pI = points[p];

        const labelList& cells = mesh_.pointCells()[pI];

        forAll(cells, c)
        {
            const label& cellN = cells[c];

            if(cellN != cellI)
            {
                if(findIndex(neighbCells, cellN) == -1)
                {
                    neighbCells.append(cells[c]);
                }
            }
        }
    }

    labelList cells;

    cells.transfer(neighbCells.shrink());

//     forAll(cells, c)
//     {
//         Info << "   cellN: " << cells[c]
//             << ", cell-centre: "
//             << mesh_.cellCentres()[cells[c]]
//             << endl;
//     }

    return cells;
}

// include cell only if it has not already been included in the DIL
template<class ParticleType>
bool Foam::cellInteractions<ParticleType>::includeCell
(
    const List<DynamicList<label> >& cells,
    const label& cellI,
    const label& cellN
)
{
    // test first if cellN has cellI
    bool include = true;

    if(findIndex(cells[cellN], cellI) != -1)
    {
        include = false;
    }

    // test also if cellI already has cellN
    if(findIndex(cells[cellI], cellN) != -1)
    {
        include = false;
    }

    return include;
}

template<class ParticleType>
bool Foam::cellInteractions<ParticleType>::interactingCells
(
    const label& cellI,
    const label& cellN,
    const scalar& offset
)
{
    bool interacting = false;

    boundedBox bb1 = cellToBoundBox(cellI);
    boundedBox bb2 = cellToBoundBox(cellN);

    bb1.inflate(offset);

    if(bb1.justOverlaps(bb2))
    {
        interacting = true;
    }

    return interacting;
}


template<class ParticleType>
Foam::boundedBox Foam::cellInteractions<ParticleType>::cellToBoundBox(const label& cellI)
{
    pointField points = boundingCellPoints(cellI);

    boundedBox bb(points, false);

//     Info << "bounding box -- min: " << bb.min() << " max: " << bb.max() << endl;

    return bb;
}

template<class ParticleType>
Foam::pointField Foam::cellInteractions<ParticleType>::boundingCellPoints
(
    const label& cellI
)
{
    const labelList& points = mesh_.cellPoints()[cellI];
    const labelList& faces = mesh_.cells()[cellI];

    label sizeOfList = points.size() + faces.size();

    pointField vectorPoints(sizeOfList, vector::zero);

    label counter = 0;

    forAll(points, p)
    {
        vectorPoints[counter] = mesh_.points()[points[p]];
        counter++;
    }

//     Info << "cell" << cellI << " points: " << vectorPoints
//          << ", faces: " << mesh_.cells()[cellI] << endl;

    forAll(faces, f)
    {
        vectorPoints[counter] = mesh_.faceCentres()[faces[f]];
        counter++;
    }

//     Info << "Face points: " << vectorFacePoints << endl;

//     Info << "points: " << vectorPoints << endl;

    return vectorPoints;
}



//- NEW
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::buildReferredCells()
{
    Info << nl << "Building referred cells..." << endl;

    //- determine source cells within interaction range of cyclic boundaries

    labelListList sourceCells(cyclics_.cyclicBoundaryModels().size());

    forAll(cyclics_.cyclicBoundaryModels(), i)
    {
        patchLayer layer
        (
            mesh_,
            rCut_,
            cyclics_.cyclicBoundaryModels()[i]->boundaryPoints(),
            -cyclics_.cyclicBoundaryModels()[i]->normal(),
            cyclics_.cyclicBoundaryModels()[i]->centroid()
        );

        const labelList& cells = layer.cells();

        // Testing
//         Info << "Cyclic boundary index = " << i << ", name = " << cyclics_.cyclicBoundaryNames()[i]
//             << endl;
//
//         forAll(cells, j)
//         {
//             Info << "cell = " << cells[j] << ", cellCentre = " << mesh_.cellCentres()[cells[j]] << endl;
//         }


        if(cells.size() > 0)
        {
            sourceCells[i].setSize(cells.size());

            forAll(cells, j)
            {
                sourceCells[i][j]=cells[j];
            }
        }

//         Info << "Test separation vector = "
//         << cyclics_.cyclicBoundaryModels()[i]->separationVector()
//         << endl;
    }

//     Info << "sourceCells = " << sourceCells << endl;

    /* 1-old
    //- each source cell is made a referred cell,
    //- referred using translation and rotational transforms
    label nRefCells = 0;

    List< DynamicList<referredCell> > referredCells(Pstream::nProcs());

    label pN = Pstream::myProcNo();

    forAll(sourceCells, i)
    {
        forAll(sourceCells[i], j)
        {
            referredCell refCell(mesh_, rCut_, sourceCells[i][j]);
            refCell.setOffsetBoundBox();

            // future include rotation
//             referredCells[i][j].setRotate
//             (
//                 rotate[patch],
//                 cyclicRotationPoints[patch],
//                 rAB[patch]
//             );

            refCell.translate(cyclics_.cyclicBoundaryModels()[i]->separationVector());
            referredCells[pN].append(refCell);
            nRefCells++;
        }
    }
    */
//     Info << "no of refCells after first cyclic boundary shift = " << nRefCells << endl;

    //- overlapping source cells

    // outer list = no of cells on mesh
    // inner list = index to cyclic boundary - i.e. size of inner list = no of boundaries the cell touches
    List<DynamicList<label> > ovSrcCells(mesh_.nCells());

    forAll(sourceCells, i)
    {
        forAll(sourceCells[i], j)
        {
            label cellI = sourceCells[i][j];

            if(findIndex(ovSrcCells[cellI], i) == -1 )
            {
                ovSrcCells[cellI].append(i);
            }

            forAll(sourceCells, k)
            {
                if(k != i)
                {
                    if(findIndex(sourceCells[k], cellI) != -1)
                    {
                        if(findIndex(ovSrcCells[cellI], k) == -1)
                        {
                            ovSrcCells[cellI].append(k);
                        }
                    }
                }
            }
        }
    }

//     Info << "ovSrcCells = " << ovSrcCells << endl;

    /*
    forAll(ovSrcCells, i)
    {
        ovSrcCells[i].shrink();

        // correction step here for 6 boundaries (i.e. single cell simulation)

        if(ovSrcCells[i].size() == 6)
        {
            DynamicList<label> mod(0);
            DynamicList<label> mod2(0);

            forAll(cyclics_.cyclicBoundaryModels(), j)
            {
                label pI = cyclics_.cyclicBoundaryModels()[j]->patchId();
                label pIN = cyclics_.cyclicBoundaryModels()[j]->patchIdN();

//                 Info << "pI = " << pI << endl;
//                 Info << "pIN = " << pIN << endl;

                if(
                    ( findIndex(mod, pI) == -1 ) &&
                    ( findIndex(mod, pIN) == -1 )
                )
                {
                    mod.append(pI);
                    mod.append(pIN);
                    mod2.append(j);
                }
            }

            mod.shrink();
            mod2.shrink();

//             Info << "mod1 = " << mod << endl;
//             Info << "mod2 = " << mod2 << endl;

            ovSrcCells[i].clear();
            ovSrcCells[i].transfer(mod2);

        }
    }
    */

//     Info << "ovSrcCells = " << ovSrcCells << endl;

    // new code for building referred cells

    label nRefCells = 0;

    List< DynamicList<referredCell> > referredCells(Pstream::nProcs());

    label pN = Pstream::myProcNo();

    forAll(ovSrcCells, i)
    {
        label nOvs=ovSrcCells[i].size();

        for (label n1 = 0; n1 < nOvs; n1++)
        {
            label k1 = ovSrcCells[i][n1];
            referredCell refCell1(mesh_, rCut_, i);
//             refCell1.setOffsetBoundBox();
            refCell1.translate(cyclics_.cyclicBoundaryModels()[k1]->separationVector());
            referredCells[pN].append(refCell1);
            nRefCells++;

            for (label n2 = 0; n2 < nOvs; n2++)
            {
                if(n2 > n1)
                {
                    label k2 = ovSrcCells[i][n2];
                    referredCell refCell2(refCell1);
                    refCell2.translate(cyclics_.cyclicBoundaryModels()[k2]->separationVector());
                    referredCells[pN].append(refCell2);
                    nRefCells++;

                    for (label n3 = 0; n3 < nOvs; n3++)
                    {
                        if(n3 > n2)
                        {
                            label k3 = ovSrcCells[i][n3];
                            referredCell refCell3(refCell2);
                            refCell3.translate(cyclics_.cyclicBoundaryModels()[k3]->separationVector());
                            referredCells[pN].append(refCell3);
                            nRefCells++;
                        }
                    }
                }
            }
        }
    }

//     Info << "number of referredCells = " << nRefCells << ", " << referredCells[pN].size() << endl;


    // check for overlapping referred cells

//     if(Pstream::parRun())
//     {
//         List< DynamicList<referredCell> > referredCells(Pstream::nProcs());
//     }
//     else
    {
//         List< DynamicList<referredCell> > referredCells(Pstream::nProcs());
//         scalar treshold = 0.01;
        scalar s = rCut_*0.2;

        List<DynamicList<referredCell> > newReferredCells(Pstream::nProcs());

        forAll(referredCells[pN], i)
        {
//             vector rI = referredCells[pN][i].midpoint();

            bool overlap = false;

            forAll(referredCells[pN], j)
            {
                if(j > i)
                {
                    referredCells[pN][j].contractII(s);

                    if(referredCells[pN][i].contains(referredCells[pN][j])) // overlap
                    {
                        overlap = true;
                    }

                    referredCells[pN][j].expandII(s);
                }
            }

            if(!overlap)
            {
                // bug 2
                // check that referred cells are not inside the domain
                referredCells[pN][i].contractII(s);

                boundedBox bb( mesh_.bounds().min(), mesh_.bounds().max());

                if(bb.contains(referredCells[pN][i])) // overlap
                {
                    referredCells[pN][i].expandII(s);
                }
                else
                {
                    referredCells[pN][i].expandII(s);
                    newReferredCells[pN].append(referredCells[pN][i]);
                }
            }
        }

//         Info << "number of new referredCells = " << newReferredCells[pN].size() << endl;

        referredCells[pN].clear();
        referredCells[pN].transfer(newReferredCells[pN]);
    }

    // now prepare referred cells with offset box
    forAll(referredCells[pN], i)
    {
        referredCells[pN][i].setOffsetBoundBox();
    }

    /* 2 old
    //- referred cells
    // for all cells on mesh
    forAll(ovSrcCells, i)
    {
        label nOvs=ovSrcCells[i].size();

        //- if nOvs = 2, => +1 corner refCell, nOvs = 3, => +4 combinations of refCells
        forAll(ovSrcCells[i], j)
        {
            forAll(ovSrcCells[i], k)
            {
                if(k > j)
                {
                    referredCell refCell(mesh_, rCut_, i);

                    refCell.setOffsetBoundBox();

                    label cyclicIndexA = ovSrcCells[i][j];

                    refCell.translate(cyclics_.cyclicBoundaryModels()[cyclicIndexA]->separationVector());

                    label cyclicIndexB = ovSrcCells[i][k];

                    refCell.translate(cyclics_.cyclicBoundaryModels()[cyclicIndexB]->separationVector());

                    referredCells[pN].append(refCell);

                    nRefCells++;
                }
            }
        }

        if(nOvs == 3)
        {
//             Info << " three overlaps !!" << endl;

            referredCell refCell(mesh_, rCut_, i);

            refCell.setOffsetBoundBox();

            forAll(ovSrcCells[i], j)
            {
                label cyclicIndex = ovSrcCells[i][j];

                refCell.translate(cyclics_.cyclicBoundaryModels()[cyclicIndex]->separationVector());
            }

            referredCells[pN].append(refCell);

            nRefCells++;
        }
    }
    */


//     Info << "no of refCells after second cyclic boundary overlap shifts = " << nRefCells << endl;

    //- find source cells within interaction of processor-only boundaries

    if(Pstream::parRun())
    {
        DynamicList<label> srcCells;

        // source cells are also internal cells
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            // just processor
            if
            (
                isA<processorPolyPatch>(patch) &&
                !isA<processorCyclicPolyPatch>(patch)
            )
            {
                nProcessorBoundaries_++;

                labelList sProcCells = getProcessorSourceCells(patch);

                forAll(sProcCells, i)
                {
                    if(findIndex(srcCells, sProcCells[i]) == -1)
                    {
                        srcCells.append(sProcCells[i]);
                    }
                }
//                 Info<< "name = " << patch.name()
//                     << ", faces = " << patch.size()
//                     <<", sourceCells = " << srcCells
//                     << endl;
            }
        }

        srcCells.shrink();

        forAll(srcCells, i)
        {
            referredCell refCell(mesh_, rCut_, srcCells[i]);

            refCell.setOffsetBoundBox();

            referredCells[pN].append(refCell);
        }
    }

    referredCells[pN].shrink();

    //- send/receive referredCells to all processors

    if (Pstream::parRun())
    {
        List<referredCell> refList(referredCells[pN].size());

        forAll(referredCells[pN], i)
        {
            refList[i] = referredCells[pN][i];
        }

        //- delete processor-referred cells on own processor (solving bug)

        referredCells[pN].clear();

        forAll(refList, i)
        {
            if(mag(refList[i].translation()) > 0)
            {
                referredCells[pN].append(refList[i]);
            }
        }

        referredCells[pN].shrink();

        //- sending
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const label proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << refList;
                }
            }
        }

        //- receiving
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<referredCell> refCellsProc;

                const label proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> refCellsProc;
                }

                forAll(refCellsProc, i)
                {
                    referredCells[p].append(refCellsProc[i]);
                }
            }
        }
    }

    //- if referredCell finds cells within boundedBox then setup interaction with processor

    DynamicList<referredCell> localRefCells;

    List<DynamicList<label> > sendSrcCells(Pstream::nProcs());
    List<DynamicList<label> > recRefIds(Pstream::nProcs());

    label count = 0;

    forAll(referredCells, p)
    {
        forAll(referredCells[p], i)
        {
            DynamicList<label> cellsN;

            forAll(mesh_.cells(), c)
            {
                if(referredCells[p][i].justOverlaps(cellToBoundBox(c)))
                {
                    cellsN.append(c);
                }
            }

            cellsN.shrink();
            List<label> neighbCells(cellsN.size());
            neighbCells.transfer(cellsN);
            referredCells[p][i].setNeighbouringCells(neighbCells);

            if(neighbCells.size() > 0)
            {
                localRefCells.append(referredCells[p][i]);
                sendSrcCells[p].append(referredCells[p][i].sourceCell());
                recRefIds[p].append(count);
                count++;
            }
        }
    }

    localRefCells.shrink();

    refCells_.setSize(localRefCells.size());

    forAll(refCells_, i)
    {
        refCells_[i] = localRefCells[i];
    }

//     Info << "refCells total = " << refCells_.size() << endl;

    // Test

//     forAll(refCells_,i)
//     {
//         if(refCells_[i].sourceCell() == 7)
//         {
//             Info << "cell centre for source cell 7," <<  mesh_.cellCentres()[7] << endl;
//             Info<< "ref cell index = " <<  i << " for source cell 7, bound box midpoint = "
//                 << refCells_[i].midpoint()
//                 << ", translation = " << refCells_[i].translation()
//                 << endl;
//
//             const labelList& realCells =refCells_[i].neighbouringCells();
//
//             forAll(realCells, j)
//             {
//                 Info << "real cell = " << realCells[j]
//                 << ", cell centre = " << mesh_.cellCentres()[realCells[j]]
//                 << endl;
//
//             }
//         }
//     }

    refCellsParticles_.setSize(refCells_.size());

    recRefIds_.setSize(Pstream::nProcs());

    forAll(recRefIds_, p)
    {
        recRefIds[p].shrink();
        recRefIds_[p].transfer(recRefIds[p]);
    }

    {
        //- sending
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            sendSrcCells[p].shrink();

            labelList srcCellsToProc(sendSrcCells[p].size());

            srcCellsToProc.transfer(sendSrcCells[p]);

            if(p != Pstream::myProcNo())
            {
                const label proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << srcCellsToProc;
                }
            }
            else // set own
            {
                sendSrcCells_[Pstream::myProcNo()].setSize(srcCellsToProc.size());
                sendSrcCells_[Pstream::myProcNo()]=srcCellsToProc;
            }
        }

        //- receiving
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                labelList srcCellsFromProc;

                const label proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> srcCellsFromProc;
                }

                sendSrcCells_[p].setSize(srcCellsFromProc.size());
                sendSrcCells_[p]=srcCellsFromProc;
            }
        }
    }

    //- build sourceCellToRefs

    List< DynamicList< label > > sourceCellToRefs(mesh_.nCells());

    forAll(refCells_, r)
    {
        const referredCell& refCellI = refCells_[r];

        const label& sourceCellI = refCellI.sourceCell();

        if(refCellI.origProcNo() == Pstream::myProcNo())
        {
            sourceCellToRefs[sourceCellI].append(r);
        }
    }

    forAll(sourceCellToRefs, r)
    {
        sourceCellToRefs[r].shrink();
        sourceCellToRefs_[r].setSize(sourceCellToRefs[r].size());
        sourceCellToRefs_[r].transfer(sourceCellToRefs[r]);
    }



    //- build inverseFRIL_

    List< DynamicList< label > > inverseFRIL(mesh_.nCells());

    forAll(refCells_, r)
    {
        const referredCell& refCellI = refCells_[r];

        const labelList& neighbours = refCellI.neighbouringCells();

        forAll(neighbours, n)
        {
            inverseFRIL[neighbours[n]].append(r);
        }
    }

    forAll(inverseFRIL, c)
    {
        inverseFRIL[c].shrink();
        inverseFRIL_[c].setSize(inverseFRIL[c].size());
        inverseFRIL_[c].transfer(inverseFRIL[c]);
    }


    Info << nl << "...done" << endl;

    // Write out referred cells
    writeReferredCells();


    // Testing

//     forAll(refCells_, i)
//     {
//         vector t = refCells_[i].midpoint();
//         refCells_[i].transformPoint(t);
//         vector m = refCells_[i].midpoint();
//
//         Info<< "rotate = " << refCells_[i].rotate()
//             << ", midpoint = " << m
//             << ", translation = " << refCells_[i].translation()
//             << ", transform midpoint = " << t
//             << endl;
//     }

}

template<class ParticleType>
Foam::labelList Foam::cellInteractions<ParticleType>::getProcessorSourceCells
(
    const polyPatch& patch
)
{
    DynamicList<label> cells;

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        scalar offset = 0.01;
        boundedBox bb = faceToBoundBox(globalFaceI, offset);

//         bb.expand(rCut_-offset);
        bb.expandII(rCut_-offset);

        forAll(mesh_.cells(), c)
        {
            if(bb.justOverlaps(cellToBoundBox(c)))
            {
                if(findIndex(cells, c) == -1)
                {
                    cells.append(c);
                }
            }
        }
    }

    cells.shrink();

    labelList sourceCells(cells.size());

    sourceCells.transfer(cells);

    return sourceCells;
}


template<class ParticleType>
Foam::boundedBox Foam::cellInteractions<ParticleType>::faceToBoundBox
(
    const label& faceI,
//     const scalar& offset,
    const scalar& smallOffset
)
{
    const labelList& nodes = mesh_.faces()[faceI];

    pointField points(nodes.size());

    forAll(nodes, p)
    {
        points[p] = mesh_.points()[nodes[p]];
    }

    vector nF = mesh_.faceAreas()[faceI];
    nF /= mag(nF);

    pointField newPoints(nodes.size()*2);

    label counter = 0;

    forAll(nodes, p)
    {
        newPoints[counter] = mesh_.points()[nodes[p]] + smallOffset*nF;
        counter++;
    }

    forAll(nodes, p)
    {
        newPoints[counter] = mesh_.points()[nodes[p]] - smallOffset*nF;
        counter++;
    }

    boundedBox bb(newPoints, false);

    return bb;
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::cellInteractions<ParticleType>::setRIPL()
{
    // clear interaction lists
    forAll(ripl_, c)
    {
        ripl_[c].clear();
    }

    forAll(inverseFRIL_, c)
    {
        const labelList& refCellIds = inverseFRIL_[c];

        forAll(refCellIds, i)
        {
            const label& r = refCellIds[i];

            forAll(refCellsParticles_[r], k)
            {
                ripl_[c].append(refCellsParticles_[r][k]);
            }
        }
    }

    forAll(ripl_, c)
    {
        ripl_[c].shrink();
    }
}




// NEW VERSION
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::setReferredParticles
(
    const List<DynamicList<ParticleType*> >& cellOccupancy
)
{
    //- clear
    referredCloud_.clear();

    forAll(refCellsParticles_, i)
    {
        refCellsParticles_[i].clear();
    }

    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    IDLList<ParticleType> ownMeshTransferList;
    labelList ownMeshParticleCount;


    //- sending
    for (label p = 0; p < Pstream::nProcs(); p++)
    {
        if(sendSrcCells_[p].size() > 0)
        {
            if(p != Pstream::myProcNo())
            {
                labelList particleCount(sendSrcCells_[p].size(), 0);
                IDLList<ParticleType> particleTransferList;

                forAll(sendSrcCells_[p], i)
                {
                    label cellI = sendSrcCells_[p][i];
                    List<ParticleType*> realParticles = cellOccupancy[cellI];
                    particleCount[i]=realParticles.size();

                    forAll(realParticles, j)
                    {
                        const ParticleType& particle = *realParticles[j];
                        particleTransferList.append(particle.clone().ptr());
                    }
                }

                UOPstream particleStream
                (
                    p,
                    pBufs
                );

                particleStream
                    << particleCount
                    << particleTransferList;
            }
            else
            {
                ownMeshParticleCount.setSize(sendSrcCells_[p].size(), 0);

                forAll(sendSrcCells_[p], i)
                {
                    label cellI = sendSrcCells_[p][i];
                    List<ParticleType*> realParticles = cellOccupancy[cellI];
                    ownMeshParticleCount[i]=realParticles.size();

                    forAll(realParticles, j)
                    {
                        const ParticleType& particle = *realParticles[j];
                        ownMeshTransferList.append(particle.clone().ptr());
                    }
                }
            }
        }
    }

    /*labelListList allNTrans(Pstream::nProcs()); // DELETED VINCENT

    pBufs.finishedSends(allNTrans);

    bool transfered;
    transfered = false;

    forAll(allNTrans, i)
    {
        forAll(allNTrans[i], j)
        {
            if (allNTrans[i][j])
            {
                transfered = true;
                break; // not sure
            }
        }
    }

//    if (!transfered)
    //{
        //break;
    //}*/


    // NEW VINCENT
    labelList allNTrans(Pstream::nProcs());
    pBufs.finishedSends(allNTrans);


    bool transfered = false;

    forAll(allNTrans, i)
    {
        if (allNTrans[i])
        {
            transfered = true;
            break;
        }
    }
    reduce(transfered, orOp<bool>());

    /*if (!transfered)
    {
        break;
    }*/
    // END NEW VINCENT


    //- receiving
    for (label p = 0; p < Pstream::nProcs(); p++)
    {
        if(recRefIds_[p].size() > 0)
        {
            if(p != Pstream::myProcNo())
            {
                //label nRec = allNTrans[p][Pstream::myProcNo()];
                label nRec = allNTrans[p]; // NEW VINCENT

                if (nRec)
                {
                    UIPstream particleStream(p, pBufs);

                    labelList particleCount(particleStream);

                    IDLList<ParticleType> newParticles
                    (
                        particleStream,
                        typename ParticleType::iNew(mesh_)
                    );

                    // NB: re coded to solve for bug of empty cells

                    label refCellIdCounter = 0;

                    typename IDLList<ParticleType>::iterator newpIter(newParticles.begin());

                    forAll(particleCount, j)
                    {
                        for(label i = 0; i < particleCount[j]; i++)
                        {
                            ParticleType& newp = newpIter();
                            label r = recRefIds_[p][refCellIdCounter];

                            refCells_[r].transformPoint(newp.position());
                            newp.transformProperties(refCells_[r].translation());
                            newp.setAsReferred();
                            referredCloud_.append(newParticles.remove(&newp));
                            refCellsParticles_[r].append(referredCloud_.last());
                            ++newpIter;
                        }

                        refCellIdCounter++;
                    }
                }
            }
            else
            {
                // NB: re coded to solve for bug of empty cells
                label refCellIdCounter = 0;

                typename IDLList<ParticleType>::iterator newpIter(ownMeshTransferList.begin());

                forAll(ownMeshParticleCount, j)
                {
                    for(label i = 0; i < ownMeshParticleCount[j]; i++)
                    {
                        ParticleType& newp = newpIter();
                        label r = recRefIds_[p][refCellIdCounter];

                        refCells_[r].transformPoint(newp.position());
                        newp.transformProperties(refCells_[r].translation());
                        newp.setAsReferred();
                        referredCloud_.append(ownMeshTransferList.remove(&newp));
                        refCellsParticles_[r].append(referredCloud_.last());
                        ++newpIter;
                    }

                    refCellIdCounter++;
                }
            }
        }
    }

    // only for debugging
//     checkForOverlaps();
}

template<class ParticleType>
void Foam::cellInteractions<ParticleType>::checkForOverlaps()
{
    DynamicList<ParticleType*> mols;

    forAllIter
    (
        typename IDLList<ParticleType>,
        referredCloud_,
        mol
    )
    {
        mols.append(&mol());
    }

    mols.shrink();

    forAll(mols, i)
    {
        vector rI = mols[i]->position();

        forAll(mols, j)
        {
            if(j>i)
            {
                vector rJ = mols[j]->position();

                if(mag(rI - rJ) < 0.01)
                {
                    FatalErrorIn("cellInteractions::checkForOverlaps Test") << nl
                        << "overlap found at position = " << rI
                        << nl << abort(FatalError);
                }
            }

        }
    }
}

//serial assumption
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::addParticle
(
    ParticleType* particle,
    DynamicList<ParticleType*>& newRefParticles,
    DynamicList<label>& refIds
)
{
    label p = Pstream::myProcNo();

    const label& cellI = particle->cell();

    labelList cS = findIndices(sendSrcCells_[p], cellI);

//     Info << "sendSrcCells[p] = "  << sendSrcCells_[p] << endl;

    if(cS.size() != 0)
    {
//         Info << "position of starting molecule = " << particle->position() << endl;

        forAll(cS, j)
        {
//             Info<< " position =" << particle->position()
//                 << ", SOURCE CELL = " << cellI << ", c = " << c
//                 << endl;
            label c = cS[j];
            label r = recRefIds_[p][c];

//             label srcCell = refCells_[r].sourceCell();
//
//             Info << "source cell check = " << srcCell << endl;

            const ParticleType& molI = *particle;
            referredCloud_.append(molI.clone().ptr());

            ParticleType* newMol = referredCloud_.last();

            refCells_[r].transformPoint(newMol->position());

            newMol->setAsReferred();

            const labelList& realCells = refCells_[r].neighbouringCells();

//             Info << "new position added at = " << newMol->position() << endl;

            forAll(realCells, rC)
            {
                ripl_[realCells[rC]].append(newMol);
                ripl_[realCells[rC]].shrink();
            }

            newRefParticles.append(newMol);
            refIds.append(r);
        }
    }

    newRefParticles.shrink();
    refIds.shrink();
}

template<class ParticleType>
void Foam::cellInteractions<ParticleType>::deleteParticle
(
    ParticleType* particle,
    DynamicList<ParticleType*>& newRefParticles,
    DynamicList<label>& refIds
)
{
    scalar tolerance = 0.01;

    label p = Pstream::myProcNo();

    const label& cellI = particle->cell();

    labelList cS = findIndices(sendSrcCells_[p], cellI);

    if(cS.size() != 0)
    {
//         Info << "position of starting molecule = " << particle->position()
//              << " cell = " << cellI
//              << endl;

        forAll(cS, j)
        {
            label c = cS[j];

            label r = recRefIds_[p][c];

            vector newPosition = particle->position();

            refCells_[r].transformPoint(newPosition);

//             Info << "source cell -> " << refCells_[r].sourceCell()
//                  << ", translate -> " << refCells_[r].translation()
//                  << ", new position = " << newPosition
//                  << endl;

            {
                typename IDLList<ParticleType>::iterator mol(referredCloud_.begin());

                for (mol = referredCloud_.begin(); mol != referredCloud_.end(); ++mol)
                {
                    ParticleType* newMol = &mol();

                    if(mag(newMol->position() - newPosition) < tolerance )
                    {
                        newRefParticles.append(newMol);
                        refIds.append(r);

//                         Info << "referred molecule found = " << newMol->position() << endl;
                    }
                }
            }

            const labelList& realCells = refCells_[r].neighbouringCells();

            forAll(realCells, rC)
            {
                DynamicList<ParticleType*> newList;

                forAll(ripl_[realCells[rC]], i)
                {
                    ParticleType* newMol = ripl_[realCells[rC]][i];

                    if(mag(newMol->position() - newPosition) < tolerance )
                    {}
                    else
                    {
                        newList.append(newMol);
                    }
                }

                ripl_[realCells[rC]].clear();
                newList.shrink();
                ripl_[realCells[rC]].transfer(newList);
            }
        }
    }

    newRefParticles.shrink();
    refIds.shrink();
}



template<class ParticleType>
void Foam::cellInteractions<ParticleType>::checkMesh(const scalar& rCut)
{
    Info << "cellInteractions: checkMesh for rCut = " << rCut << endl;

    const polyMesh& mesh = mesh_;

    scalar tolerance = 0.001*rCut;

    for (label c = 0; c < mesh.nCells(); c++)
    {
        // by bounding box

        boundedBox bb(boundingCellPoints(c), false);

        if(bb.span().x() < (rCut - tolerance))
        {
            FatalErrorIn("cellInteractions::checkMesh") << nl
                << "WARNING: Test for bound box failed! " << nl
                << "Cell: " << c << " cell-centre: " << mesh.cellCentres()[c]
                << " distance of x: " << bb.span().x()
                << nl << " Modify mesh to make cells greater or equal to rCut = " << rCut
                << nl << " OR else use an unoptimised mesh configuration, by including:"
                << nl <<" optimisedMesh       no; " << nl << "in the system/potentialsDict."
                << nl << abort(FatalError);

        }

        if(bb.span().y() < (rCut - tolerance))
        {
            FatalErrorIn("cellInteractions::checkMesh") << nl
                << "WARNING: Test for bound box failed! " << nl
                << "Cell: " << c << " cell-centre: " << mesh.cellCentres()[c]
                << " distance of y: " << bb.span().y()
                << nl << " Modify mesh to make cells greater or equal to rCut = " << rCut
                << nl << " OR else use an unoptimised mesh configuration, by including:"
                << nl <<" optimisedMesh       no; " << nl << "in the system/potentialsDict."
                << nl << abort(FatalError);
        }

        if(bb.span().z() < (rCut - tolerance))
        {
            FatalErrorIn("cellInteractions::checkMesh") << nl
                << "WARNING: Test for bound box failed! " << nl
                << "Cell: " << c << " cell-centre: " << mesh.cellCentres()[c]
                << " distance of z: " << bb.span().z()
                << nl << " Modify mesh to make cells greater or equal to rCut = " << rCut
                << nl << " OR else use an unoptimised mesh configuration, by including:"
                << nl <<" optimisedMesh       no; " << nl << "in the system/potentialsDict."
                << nl << abort(FatalError);
        }

        // by volume
        scalar V = mesh.cellVolumes()[c];

        scalar r = Foam::pow(V, (1.0/3.0) );

        if(r < (rCut - tolerance))
        {
            FatalErrorIn("cellInteractions::checkMesh") << nl
                << "WARNING: Test for volume failed! " << nl
                << "Cell: " << c << " cell-centre: " << mesh.cellCentres()[c]
                << " avarage r = " << r
                << nl << " Modify mesh to make cells greater or equal to rCut = " << rCut
                << nl << " OR else use an unoptimised mesh configuration, by including:"
                << nl <<" optimisedMesh       no; " << nl << "in the system/potentialsDict."
                << nl << abort(FatalError);
        }
    }
}






// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class ParticleType>
bool Foam::cellInteractions<ParticleType>::write()
{
    return write_;
}

//- visualisation of referred cells using VMD - only for future debugging
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::writeReferredCells()
{
    if(writeReferredCells_)
    {
        if(Pstream::parRun())
        {
            Info << "Writing out referred cells" << endl;

            std::string s;
            std::stringstream out;
            out << Pstream::myProcNo();
            s = out.str();

            {
                fileName fName1("referredCells_procNo_"+s+"_cyclic.xmol");
                fileName fName2("referredCells_procNo_"+s+"_proc.xmol");
                fileName fName3("referredCells_procNo_"+s+"_cyclic_RU.xmol");
                fileName fName4("referredCells_procNo_"+s+"_proc_RU.xmol");

                label nProcCells = 0;
                label nCycCells = 0;

                forAll(refCells_, i)
                {
                    if(mag(refCells_[i].translation()) > SMALL)
                    {
                        nCycCells++;
                    }
                    else
                    {
                        nProcCells++;
                    }
                }

                OFstream os1(fName1);
                OFstream os2(fName2);
                OFstream os3(fName3);
                OFstream os4(fName4);

                os1 << nCycCells << nl << "cyclic referredCells positions in angstroms" << nl;
                os2 << nProcCells << nl << "processor referredCells positions in angstroms" << nl;
                os3 << nCycCells << nl << "cyclic referredCells positions in reduced units" << nl;
                os4 << nProcCells << nl << "processor referredCells positions in reduced units" << nl;

                forAll(refCells_, i)
                {
                    vector vA = refCells_[i].midpoint()*redUnits_.refLength()*1e10;
                    vector vRU = refCells_[i].midpoint();

                    if(mag(refCells_[i].translation()) > SMALL)
                    {
                        os1 << "C "
                            << ' ' << vA.x()
                            << ' ' << vA.y()
                            << ' ' << vA.z()
                            << nl;

                        os3 << "C "
                            << ' ' << vRU.x()
                            << ' ' << vRU.y()
                            << ' ' << vRU.z()
                            << nl;
                    }
                    else
                    {
                        os2 << "P "
                            << ' ' << vA.x()
                            << ' ' << vA.y()
                            << ' ' << vA.z()
                            << nl;

                        os4 << "P "
                            << ' ' << vRU.x()
                            << ' ' << vRU.y()
                            << ' ' << vRU.z()
                            << nl;
                    }
                }
            }
        }
        else
        {
            Info << "Writing out referred cells" << endl;

            fileName fName1("referredCells_cyclic.xmol");
            fileName fName2("referredCells_cyclic_RU.xmol");

            OFstream os1(fName1);
            OFstream os2(fName2);

            os1 << refCells_.size() << nl << "cyclic referredCells positions in angstroms" << nl;
            os2 << refCells_.size() << nl << "cyclic referredCells positions in RU" << nl;

            forAll(refCells_, i)
            {
                vector v = refCells_[i].midpoint()*redUnits_.refLength()*1e10;

                os1 << "C "
                    << ' ' << v.x()
                    << ' ' << v.y()
                    << ' ' << v.z()
                    << nl;

                vector v2 = refCells_[i].midpoint();

                os2 << "C "
                    << ' ' << v2.x()
                    << ' ' << v2.y()
                    << ' ' << v2.z()
                    << nl;
            }
        }
    }
}



/*
template<class ParticleType>
void Foam::cellInteractions<ParticleType>::writeReferredCloud()
{
    if(writeCloud_)
    {
            const Time& runTime = mesh_.time();
//
//         if(runTime.outputTime())
//         {
            Info << "Writing out referred cloud" << endl;

            fileName timePath(runTime.path()/runTime.timeName()/"lagrangian"/fieldName_+"_ReferredCloud");

            if (!isDir(timePath))
            {
                mkDir(timePath);
            }

            fileName fName(timePath/"positions");
            fileName idName(timePath/"id");

            OFstream str(fName);
            OFstream str2(idName);

            str << " FoamFile" << nl << "{" << nl
                <<"     version     2.0;" << nl
                <<"     format      ascii;" << nl
                <<"     class       Cloud<molecule>;" << nl
                <<"     location    \"0.1\";" << nl
                <<"     object      positions;" << nl
                << nl << "}" << nl << endl;

            str << referredCloud_.size() << nl << "(" << nl;

            str2 << " FoamFile" << nl << "{" << nl
                <<"     version     2.0;" << nl
                <<"     format      ascii;" << nl
                <<"     class       labelField" << nl
                <<"     location    \"0.1\";" << nl
                <<"     object      id;" << nl
                << nl << "}" << nl << endl;

            str2 << referredCloud_.size() << nl << "(" << nl;

            forAllIter
            (
                typename IDLList<ParticleType>,
                referredCloud_,
                mol
            )
            {
                str << mol().position() << " " << -1 << nl;
                str2 << mol().id() << nl;
            }

            str << ")" << nl;
            str2 << ")" << nl;
//         }
    }
}*/



// ************************************************************************* //
