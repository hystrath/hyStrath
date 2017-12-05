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

#include "matteoliAndMansooriMod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(matteoliAndMansooriMod, 0);

addToRunTimeSelectionTable(rdfModel, matteoliAndMansooriMod, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
matteoliAndMansooriMod::matteoliAndMansooriMod
(
//     Time& t,
    const dictionary& dict
)
:
    rdfModel(/*t,*/ dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    name_(propsDict_.lookup("distributionName")),
    T_(readScalar(propsDict_.lookup("T"))),
    density_(readScalar(propsDict_.lookup("density"))),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    rMax_(readScalar(propsDict_.lookup("rMax"))),
    binWidth_(rMax_/scalar(nBins_-1)),
    g_(nBins_, 0.0),
    r_(nBins_, 0.0)
/*
    h_(0.0),
    m_(0.0),
    gd_(0.0),
    lambda_(0.0),
    alpha_(0.0),
    beta_(0.0),
    theta_(0.0)*/
{
    //- set r_

    for(label i = 0; i < nBins_; i++)
    {
       r_[i] = (0.5 + scalar(i)) * binWidth_;
    }


    //- set g_

//     const dictionary& pieceWiseDict(propsDict_.subDict("piecewiseFunction"));

    PtrList<entry> functionList(propsDict_.lookup("functions"));

    forAll(functionList, f)
    {
        const entry& function = functionList[f];
    
        const dictionary& functionDict = function.dict();

        bool tailPartTrue = false;

        tailPartTrue = Switch(functionDict.lookup("tailPart"));

        const scalar rS = readScalar(functionDict.lookup("startPoint"));        
        const scalar rE = readScalar(functionDict.lookup("endPoint"));        

        const dictionary& coeffdict(functionDict.subDict("coefficients"));
     
        if(tailPartTrue)
        {
            scalar h = readScalar(coeffdict.lookup("h"));
            scalar m = readScalar(coeffdict.lookup("m"));
            scalar gd = readScalar(coeffdict.lookup("gd"));
            scalar lambda = readScalar(coeffdict.lookup("lambda"));
            scalar alpha = readScalar(coeffdict.lookup("alpha"));
            scalar beta = readScalar(coeffdict.lookup("beta"));

            for(label i = 0; i < g_.size(); i++)
            {
                if
                (
                    (rS <= r_[i]) &&
                    (r_[i] < rE)
                )
                {
                    g_[i] = tailPart(r_[i]/h,m,gd,lambda,alpha,beta);
                }
            }
        }
        else
        {
            scalar h = readScalar(coeffdict.lookup("h"));
            scalar gd = readScalar(coeffdict.lookup("gd"));
            scalar theta = readScalar(coeffdict.lookup("theta"));

            for(label i = 0; i < g_.size(); i++)
            {
                if
                (
                    (rS <= r_[i]) &&
                    (r_[i] < rE)
                )
                {
                    g_[i] = initialPart(r_[i]/h,gd,theta);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

matteoliAndMansooriMod::~matteoliAndMansooriMod()
{}


scalar matteoliAndMansooriMod::tailPart
(
    const scalar& y,
    const scalar& m,
    const scalar& gd,
    const scalar& lambda,
    const scalar& alpha,
    const scalar& beta
)
{
    return 1.0 + pow(y, -m)*(gd-1.0-lambda) + ((y - 1.0+lambda)*exp(-alpha*(y-1.0))*cos(beta*(y-1.0)) )/y;
}


scalar matteoliAndMansooriMod::initialPart
(

    const scalar& y,
    const scalar& gd,
    const scalar& theta
)
{
    return gd*exp(-theta*sqr(y-1));
}

/*
const scalarField& matteoliAndMansooriMod::g() const
{
    return g_;
}


const scalarField& matteoliAndMansooriMod::r() const
{
    return r_;
}

const scalar& matteoliAndMansooriMod::binWidth() const
{
    return binWidth_;
}*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void matteoliAndMansooriMod::setRDF(radialDistribution& rdf, const Time& runTime)
{
    rdf.setRdf
    (
        name_,
        g_,
        r_
    );
}


} // End namespace Foam

// ************************************************************************* //
