/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "potential.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::potential::setSiteIdList(const IOdictionary& moleculePropertiesDict)
{
    DynamicList<word> siteIdList;

    DynamicList<word> pairPotentialSiteIdList;

    dictionary moleculeProperties
    (
        moleculePropertiesDict.subDict("moleculeProperties")
    );

    forAll(idList_, i)
    {
        const word& id(idList_[i]);

        if (!moleculeProperties.found(id))
        {
            FatalErrorIn("potential::setSiteIdList(const dictionary&)")
                << id << " molecule subDict not found"
                << nl << abort(FatalError);
        }

        const dictionary& molDict(moleculeProperties.subDict(id));

        List<word> siteIdNames = molDict.lookup("siteIds");

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            if(findIndex(siteIdList, siteId) == -1)
            {
                siteIdList.append(siteId);
            }
        }

        List<word> pairPotSiteIds = molDict.lookup("pairPotentialSiteIds");

        forAll(pairPotSiteIds, sI)
        {
            const word& siteId = pairPotSiteIds[sI];

            if (findIndex(siteIdNames, siteId) == -1)
            {
                FatalErrorIn("potential::setSiteIdList(const dictionary&)")
                    << siteId << " in pairPotentialSiteIds is not in siteIds: "
                    << siteIdNames << nl << abort(FatalError);
            }

            if (findIndex(pairPotentialSiteIdList, siteId) == -1)
            {
                pairPotentialSiteIdList.append(siteId);
            }
        }
    }

    nPairPotIds_ = pairPotentialSiteIdList.size();

    forAll(siteIdList, aSIN)
    {
        const word& siteId = siteIdList[aSIN];

        if (findIndex(pairPotentialSiteIdList, siteId) == -1)
        {
            pairPotentialSiteIdList.append(siteId);
        }
    }

    siteIdList_.transfer(pairPotentialSiteIdList);
}


void Foam::potential::readPotentialDict()
{
    Info<< nl <<  "Reading potential dictionary:" << endl;

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    setIdList(moleculePropertiesDict);

    setSiteIdList(moleculePropertiesDict);

    List<word> pairPotentialSiteIdList
    (
        SubList<word>(siteIdList_, nPairPotIds_)
    );

    Info<< nl << "Unique site ids found: " << siteIdList_
        << nl << "Site Ids requiring a pair potential: "
        << pairPotentialSiteIdList
        << endl;

    IOdictionary potentialDict
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    potentialEnergyLimit_ = readScalar
    (
        potentialDict.lookup("potentialEnergyLimit")
    );


    if(redUnits_.runReducedUnits())
    {
        potentialEnergyLimit_ /= redUnits_.refEnergy();
    }

     //************ modified the way removal ordering is done *****/
    if (potentialDict.found("removalOrder"))
    {
        List<word> remOrd = potentialDict.lookup("removalOrder");

        DynamicList<label> removalOrder(0);

        forAll(remOrd, rO)
        {
            const label id = findIndex(idList_, remOrd[rO]);

            if( (id != -1) && (findIndex(removalOrder, id) == -1) )
            {
                removalOrder.append(id);
            }
        }

        forAll(idList_, i)
        {
            if(findIndex(removalOrder, i) == -1)
            {
                removalOrder.append(i);
            }
        }

        removalOrder_.transfer(removalOrder);
    }

    // *************************************************************************
    // Pair potentials

    if (!potentialDict.found("pair"))
    {
        FatalErrorIn("potential::readPotentialDict()")
            << "pair potential specification subDict not found"
            << abort(FatalError);
    }

    const dictionary& pairDict = potentialDict.subDict("pair");

    pairPotentials_.buildPotentials
    (
        pairPotentialSiteIdList,
        pairDict,
        mesh_,
        redUnits_
    );

    // *************************************************************************
    // Tether potentials

    if (tetherIdList_.size())
    {
        if (!potentialDict.found("tether"))
        {
            FatalErrorIn("potential::readPotentialDict()")
                << "tether potential specification subDict not found"
                << abort(FatalError);
        }

        const dictionary& tetherDict = potentialDict.subDict("tether");

        tetherPotentials_.buildPotentials
        (
            siteIdList_,
            tetherDict,
            tetherIdList_,
            redUnits_
        );
    }

    // *************************************************************************
	// External Forces

	gravity_ = vector::zero;

	if (potentialDict.found("external"))
	{
		Info<< nl << "Reading external forces:" << endl;

		const dictionary& externalDict = potentialDict.subDict("external");

		// gravity
		externalDict.readIfPresent("gravity", gravity_);

		Info<< nl << tab << "gravity = " << gravity_ << endl;
	}
}

void Foam::potential::setIdList
(
    const IOdictionary& moleculePropertiesDict
)
{
    // read in id-list from moleculePropertiesDict

    const List<word> molecules (moleculePropertiesDict.lookup("idList"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
        else
        {
            Info << "WARNING: YOU HAVE DEFINED MORE THAN ONCE MOLECULE: "
                 << moleculeName << nl << " in moleculeProperties dict."
                 << endl;
        }
    }

    idList_.transfer(moleculesReduced);

    // tethered molecules

    const List<word> tetheredMolecules (moleculePropertiesDict.lookup("tetherIdList"));

    DynamicList<word> tetheredMoleculesReduced(0);

    forAll(tetheredMolecules, i)
    {
        const word& moleculeName(tetheredMolecules[i]);

        if(findIndex(tetheredMoleculesReduced, moleculeName) == -1)
        {
            if(findIndex(idList_, moleculeName) != -1)
            {
                tetheredMoleculesReduced.append(moleculeName);
            }
            else
            {
                FatalErrorIn("potential::setIdList()")
                    << "Tethered molecule " << moleculeName
                    << " not specified in moleculePropertiesDict/idList"
                    << abort(FatalError);
            }
        }
    }

    tetherIdList_.transfer(tetheredMoleculesReduced);

    setSiteIdList(moleculePropertiesDict);
}


void Foam::potential::potential::readMdInitialiseDict
(
    const IOdictionary& mdInitialiseDict,
    IOdictionary& idListDict
)
{
    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    DynamicList<word> idList;

    DynamicList<word> tetherSiteIdList;

    forAll(mdInitialiseDict.toc(), zone)
    {
        const dictionary& zoneDict = mdInitialiseDict.subDict
        (
            mdInitialiseDict.toc()[zone]
        );

        List<word> latticeIds
        (
            zoneDict.lookup("latticeIds")
        );

        forAll(latticeIds, i)
        {
            const word& id = latticeIds[i];

            if (!moleculePropertiesDict.found(id))
            {
                FatalErrorIn
                (
                    "potential::readMdInitialiseDict"
                    "("
                        "const IOdictionary&, "
                        "IOdictionary&"
                    ")"
                )   << "Molecule type " << id
                    << " not found in moleculeProperties dictionary." << nl
                    << abort(FatalError);
            }

            if (findIndex(idList,id) == -1)
            {
                idList.append(id);
            }
        }

        List<word> tetherSiteIds
        (
            zoneDict.lookup("tetherSiteIds")
        );

        forAll(tetherSiteIds, t)
        {
            const word& tetherSiteId = tetherSiteIds[t];

            bool idFound = false;

            forAll(latticeIds, i)
            {
                if (idFound)
                {
                    break;
                }

                const word& id = latticeIds[i];

                List<word> siteIds
                (
                    moleculePropertiesDict.subDict(id).lookup("siteIds")
                );

                if (findIndex(siteIds, tetherSiteId) != -1)
                {
                    idFound = true;
                }
            }

            if (idFound)
            {
                tetherSiteIdList.append(tetherSiteId);
            }
            else
            {
                FatalErrorIn
                (
                    "potential::readMdInitialiseDict"
                    "("
                        "const IOdictionary&, "
                        "IOdictionary&"
                    ")"
                )   << "Tether id  " << tetherSiteId
                    << " not found as a site of any molecule in zone." << nl
                    << abort(FatalError);
            }
        }
    }

    idList_.transfer(idList);

    tetherSiteIdList.shrink();

    idListDict.add("idList", idList_);

    idListDict.add("tetherSiteIdList", tetherSiteIdList);

    setSiteIdList(moleculePropertiesDict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct to run MD simulation (idList is read in from constant dir)
Foam::potential::potential
(
    const polyMesh& mesh,
    const reducedUnits& rU
)
:
    mesh_(mesh),
    redUnits_(rU)
{
    readPotentialDict();
}

//- Construct for mdInitialise (idList built from all molecules defined in moleculesProperties)
Foam::potential::potential
(
    const polyMesh& mesh,
    const reducedUnits& rU,
    const IOdictionary& moleculePropertiesDict
)
:
    mesh_(mesh),
    redUnits_(rU)
{
    readPotentialDict();
}

Foam::potential::potential
(
    const polyMesh& mesh,
    const IOdictionary& mdInitialiseDict,
    IOdictionary& idListDict
)
:
    mesh_(mesh),
	redUnits_(reducedUnits())
{
	readMdInitialiseDict(mdInitialiseDict, idListDict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::potential::~potential()
{}


// ************************************************************************* //
