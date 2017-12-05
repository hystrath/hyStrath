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
    distributePoints

Description

\*----------------------------------------------------------------------------*/

#include "distributePoints.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distributePoints::distributePoints
(
    cachedRandomMD& rndGen
)
:
    boundedBox(),
    rndGen_(rndGen)    
{
    
    
}

distributePoints::distributePoints
(
    const boundedBox& bb,
    cachedRandomMD& rndGen
)
:
    boundedBox(bb),
    rndGen_(rndGen)
{}
    
// // Null constructor
// distributePoints::distributePoints
// (
// 	const label rndSeed,
// 	const label rndCacheSize
// )
// :
//     boundedBox(),
//     rndGen_(rndSeed, rndCacheSize) //Initialise cachedRandomMD object with provided size and seed
// {
//     
//     
// }
// 
// distributePoints::distributePoints
// (
//     const boundedBox& bb,
// 	const label rndSeed,
// 	const label rndCacheSize
// )
// :
//     boundedBox(bb),
//     rndGen_(rndSeed, rndCacheSize) //Initialise cachedRandomMD object with provided size and seed
// {}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

distributePoints::~distributePoints()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector distributePoints::randomPoint()
{
    vector pointI
    (
        this->span().x()*rndGen_.sample01<scalar>(),
        this->span().y()*rndGen_.sample01<scalar>(),
        this->span().z()*rndGen_.sample01<scalar>()
    );
    
    pointI += this->min();

    return pointI;
}

void distributePoints::setBoundedBox(const boundedBox& bb)
{
    this->min() = bb.min();
    this->max() = bb.max();
}


List<vector> distributePoints::uniform
(
    const label& nPoints
)
{
    vectorField sites(nPoints, vector::zero);
        
    if(nPoints == 1)
    {
        sites[0] = this->midpoint();
    }
    else if(nPoints > 1)
    {
        scalar Lx = this->span().x();
        scalar Ly = this->span().y();
        scalar Lz = this->span().z();
        
        scalar vol = this->volume();
        scalar density = nPoints/vol;
        scalar spacing = Foam::pow((1.0/density), (1.0/3.0));
        
        label nPtsX = label((Lx/spacing) + 0.5);
        label nPtsY = label((Ly/spacing) + 0.5);
        label nPtsZ = label((Lz/spacing) + 0.5);        
        
        label nPointsNew = nPtsX*nPtsY*nPtsZ;
        
        scalar A = 0.0;
        
        while(nPointsNew < nPoints)
        {
            A += 1.0;
            
            nPtsX = label((Lx/spacing) + A);
            nPtsY = label((Ly/spacing) + A);
            nPtsZ = label((Lz/spacing) + A);
            nPointsNew = nPtsX*nPtsY*nPtsZ;
        }
        
        scalar dRX = Lx/nPtsX;
        scalar dRY = Ly/nPtsY;
        scalar dRZ = Lz/nPtsZ;
        
        vectorField newSites(nPointsNew, vector::zero);
        
        label c = 0;
        
        for (label i = 0; i < nPtsX; i++)
        {
            for (label j = 0; j < nPtsY; j++)
            {
                for (label k = 0; k < nPtsZ; k++)
                {
                    vector x = vector(1, 0, 0)*(i*dRX + 0.5*dRX);
                    vector y = vector(0, 1, 0)*(j*dRY + 0.5*dRY);
                    vector z = vector(0, 0, 1)*(k*dRZ + 0.5*dRZ);
                    
                    newSites[c] = x + y + z + this->min();
                    c++;
                }
            }
        }
        
        // NEW RANDOM METHOD FOR SELECTING POINTS FROM POSSIBLE LIST
        DynamicList<vector> options;
        
        forAll(newSites, i)
        {
            options.append(newSites[i]);
            
        }
        
//         Pout << "distribute points = " << options.size()-1 << endl;
        
        forAll(sites, i)
        {
            if(options.size() > 1)
            {
                label n = rndGen_.integer(0, (options.size()-1));
                
                sites[i]= options[n];
                
                List<vector> transfOptions(options.size() -1);
                
                label c = 0;
                
                forAll(options, j)
                {
                    if(j != n)
                    {
                        transfOptions[c]=options[j];
                        c++;
                    }
                }
                
                options.clear();
                
                forAll(transfOptions, k)
                {
                    options.append(transfOptions[k]);
                }
                
                //options.shrink();
            }
            else if (options.size() == 1)
            {
                sites[i] = options[0];
            }
            else
            {
                Info << "WARNING - ERROR" << endl;
            }
        }
    }

    return sites;    
}


void distributePoints::testUniform()
{
    label noOfPoints= 1000;
    
    List<scalar> rMagP(noOfPoints, 0.0);
    
    for (label j = 1; j < noOfPoints; j++)
    {
        List<vector> r = uniform(j);
        
        scalar rMagAv = 0.0;
        
        List<scalar> rMags(r.size(), 0.0);
        
        forAll(r, i)
        {
            const vector& rI = r[i];
            
            scalar rMag = GREAT;
            
            forAll(r, j)
            {
                if(j != i)
                {
                    const vector& rJ = r[j];
                    
                    scalar rIJMag = Foam::mag(rI - rJ);
                    
                    if(rIJMag < rMag)
                    {
                        rMag = rIJMag;
                    }
                }
            }
            
            // also do boundaries
            
            rMagAv += rMag;
            
            rMags[i]=rMag;
        }
        
        rMagAv /= j;
        
        rMagP[j] = rMagAv;
        
        Info << j << " " << rMagAv << endl;
        
        {
            std::string s;
            std::stringstream out;
            out << j;
            s = out.str();              
            
            fileName fName2("distribution_"+s+"_RU.xmol");  //ParaFOAM

            OFstream os2(fName2);
            
            os2 << r.size() << nl << "site positions in reduced units" << nl;    
            
            forAll(r, i)
            {
                os2 << "X"
                        << ' ' << r[i].x()
                        << ' ' << r[i].y()
                        << ' ' << r[i].z()
                        << nl;
            }
        }
    }
}

void distributePoints::outputProperties(const List<vector>& r)
{
    List<scalar> rMags(r.size(), 0.0);
    
    forAll(r, i)
    {
        const vector& rI = r[i];
        
        scalar rMag = GREAT;
        
        forAll(r, j)
        {
            if(j != i)
            {
                const vector& rJ = r[j];
                
                scalar rIJMag = Foam::mag(rI - rJ);
                
                if(rIJMag < rMag)
                {
                    rMag = rIJMag;
                }
            }
        }
        
        // also do boundaries
        rMags[i]=rMag;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
