/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    selectIdPairs

Description

\*----------------------------------------------------------------------------*/

#include "selectIdPairs.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
selectIdPairs::selectIdPairs()
:
    epsilon_(),
    sigma_(),
    alpha_(),
    rho_(),
    C_(),
    lennardJones_(false),
    buckinghamPotential_(false),
    nIds_(0)
{}



selectIdPairs::selectIdPairs
(
    const polyMesh& mesh,
    const potential& pot
)
:
    epsilon_(),
    sigma_(),
    alpha_(),
    rho_(),
    C_(),
    lennardJones_(false),
    buckinghamPotential_(false),
    nIds_(0)
{
    setConfiguration(mesh, pot);
}




void selectIdPairs::setConfiguration
(
    const polyMesh& mesh,
    const potential& pot
)
{
    IOdictionary potentialDict
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

    //====================================================================
    // information about pairwise interactions given in the potentialDict:
    //====================================================================
    const dictionary& pairDict = potentialDict.subDict("pair");

    
    //==================================================
    // Display the site IDs found in the potentialDict:
    //==================================================
    Info << "OpenFOAM Site IDs:" << endl;
    Info << pot.siteIdList() << endl;


    //========================================
    // Build the potential parameter matrices:
    //========================================
    //const List<word>& idList = pot.idList(); removed line********************************
    //changed line*************************************************************************
    nIds_ = pot.siteIdList().size();

    epsilon_.setSize(nIds_);
    sigma_.setSize(nIds_);

    for (label i = 0; i < nIds_; ++i)
    {
        epsilon_[i].setSize(nIds_, 0.0);
        sigma_[i].setSize(nIds_, 0.0);
    }

    for (label a = 0; a < nIds_; ++a)
    {
        word idA = pot.siteIdList()[a];

        for (label b = a; b < nIds_; ++b)
        {
            word idB = pot.siteIdList()[b];

            word pairPotentialName;

            if (a == b)
            {
                if (pairDict.found(idA + "-" + idB))
                {
                    pairPotentialName = idA + "-" + idB;
                }
                else
                {
                    /*FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " not found"
                        << nl << abort(FatalError);*/
                    pairPotentialName = "None";  // potential pair is not defined in the potentialDict
                }
            }
            else
            {
                if (pairDict.found(idA + "-" + idB))
                {
                    pairPotentialName = idA + "-" + idB;
                }

                else if (pairDict.found(idB + "-" + idA))
                {
                    pairPotentialName = idB + "-" + idA;
                }

                else
                {
                    /*FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " or "
                        << idB << "-" << idA << " not found"
                        << nl << abort(FatalError);*/
                    pairPotentialName = "None";  // potential pair is not defined in the potentialDict
                }

                if
                (
                    pairDict.found(idA+"-"+idB)
                 && pairDict.found(idB+"-"+idA)
                )
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " and "
                        << idB << "-" << idA << " found multiple definition"
                        << nl << abort(FatalError);
                }
            }

            if(pairPotentialName=="None")
            {
               Info << idA << "-" << idB <<  ": Interaction type - No interaction." << endl;
               epsilon_[a][b] = 0;
               epsilon_[b][a] = 0;
               sigma_[a][b]   = 0;
               sigma_[b][a]   = 0;
            }
            else
            {

                const dictionary& pairPotentialDict = pairDict.subDict(pairPotentialName);

                 word pairPotentialTypeName(pairPotentialDict.lookup("pairPotential"));
 
                if(pairPotentialTypeName == "noInteraction")
                {
                    Info << idA << "-" << idB <<  ": Interaction type - No interaction." << endl;
                    epsilon_[a][b] = 0;
                    epsilon_[b][a] = 0;
                    sigma_[a][b]   = 0;
                    sigma_[b][a]   = 0;
                }
                else
                {
                    Info << idA << "-" << idB << ": Interaction type - " << pairPotentialTypeName << endl;
                    const dictionary& lJDict = pairPotentialDict.subDict("lennardJonesCoeffs");

                    scalar sigma = readScalar(lJDict.lookup("sigma"));
                    scalar epsilon = readScalar(lJDict.lookup("epsilon"));

                    epsilon_[a][b] = epsilon;
                    epsilon_[b][a] = epsilon;
                    sigma_[a][b] = sigma;
                    sigma_[b][a] = sigma;
                }
            }
//             pairMatrix_[a][b] = 
//             pairMatrix_[b][a] = 
        }
    }

    Info << "===================================" << endl;
    Info << "Select ID pairs (pair potentials): " << endl;
    Info << "===================================" << endl;
    Info << "epsilon: " << epsilon_ << endl;
    Info << "sigma: " << sigma_ << endl;
    Info << "===================================" << endl << endl;
}

/*
void selectIdPairs::setConfiguration
(
    const polyMesh& mesh,
    const potential& pot
)
{
    IOdictionary potentialDict
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

    const dictionary& pairDict = potentialDict.subDict("pair");

    const List<word>& idList = pot.siteIdList();

    nIds_ = pot.siteIdList().size();

    // lennard jones parameters
    epsilon_.setSize(nIds_);
    sigma_.setSize(nIds_);

    
    // buckingham potential parameters
    alpha_.setSize(nIds_);
    rho_.setSize(nIds_);
    C_.setSize(nIds_);    

    
    
    for (label i = 0; i < nIds_; ++i)
    {
        epsilon_[i].setSize(nIds_, 0.0);
        sigma_[i].setSize(nIds_, 0.0);

        alpha_[i].setSize(nIds_, 0.0);
        rho_[i].setSize(nIds_, 0.0);
        C_[i].setSize(nIds_, 0.0);
    }

    for (label a = 0; a < nIds_; ++a)
    {
        word idA = idList[a];

        for (label b = a; b < nIds_; ++b)
        {
            word idB = idList[b];

            word pairPotentialName;

            if (a == b)
            {
                if (pairDict.found(idA + "-" + idB))
                {
                    pairPotentialName = idA + "-" + idB;
                }
                else
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " not found"
                        << nl << abort(FatalError);
                }
            }
            else
            {
                if (pairDict.found(idA + "-" + idB))
                {
                    pairPotentialName = idA + "-" + idB;
                }

                else if (pairDict.found(idB + "-" + idA))
                {
                    pairPotentialName = idB + "-" + idA;
                }

                else
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " or "
                        << idB << "-" << idA << " not found"
                        << nl << abort(FatalError);
                }

                if
                (
                    pairDict.found(idA+"-"+idB)
                 && pairDict.found(idB+"-"+idA)
                )
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " and "
                        << idB << "-" << idA << " found multiple definition"
                        << nl << abort(FatalError);
                }
            }

            const dictionary& pairPotentialDict = pairDict.subDict(pairPotentialName);
            
            const word potentialName = pairPotentialDict.lookup("pairPotential");
            
            if(potentialName == "lennardJones")
            {
                lennardJones_ = true;
                
                const dictionary& lJDict = pairPotentialDict.subDict("lennardJonesCoeffs");

                scalar sigma = readScalar(lJDict.lookup("sigma"));
                scalar epsilon = readScalar(lJDict.lookup("epsilon"));

                epsilon_[a][b] = epsilon;
                epsilon_[b][a] = epsilon;
                sigma_[a][b] = sigma;
                sigma_[b][a] = sigma;
            }
            else if (potentialName == "buckinghamPotential")
            {
                buckinghamPotential_ = true;
                
                const dictionary& lJDict = pairPotentialDict.subDict("buckinghamPotentialCoeffs");

                scalar alpha = readScalar(lJDict.lookup("alpha"));
                scalar rho = readScalar(lJDict.lookup("rho"));
                scalar C = readScalar(lJDict.lookup("C"));                

                alpha_[a][b] = alpha;
                alpha_[b][a] = alpha;

                rho_[a][b] = rho;
                rho_[b][a] = rho;

                C_[a][b] = C;
                C_[b][a] = C;
            }
        }
    }

    if(lennardJones_)
    {
        Info<< "Lennard Jones parameters : " 
            << nl << "epsilon: " << epsilon_
            << nl << "sigma: " << sigma_
            << endl;
    }
 
    if(buckinghamPotential_)
    {
        Info<< "Buckingham potential parameters : " 
            << nl << "alpha: " << alpha_
            << nl << "rho: " << rho_
            << nl << "C: " << C_
            << endl; 
    }
}
*/


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

selectIdPairs::~selectIdPairs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
