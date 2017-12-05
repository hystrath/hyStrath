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

#include "matteoliAndMansoori.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(matteoliAndMansoori, 0);

addToRunTimeSelectionTable(rdfModel, matteoliAndMansoori, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
matteoliAndMansoori::matteoliAndMansoori
(
//     Time& t,
    const dictionary& dict
)
:
    rdfModel(/*t,*/ dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    name_(propsDict_.lookup("distributionName")),
    T_(readScalar(propsDict_.lookup("temperature"))),
    density_(readScalar(propsDict_.lookup("density"))),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    rMax_(readScalar(propsDict_.lookup("rMax"))),
    binWidth_(rMax_/(scalar(nBins_))),
    g_(nBins_, 0.0),
    r_(nBins_, 0.0),
    reducedUnits_(true),
    readCoefficients_(false),
    h_(0.0),
    m_(0.0),
    gd_(0.0),
    lambda_(0.0),
    alpha_(0.0),
    beta_(0.0),
    theta_(0.0)
{
    //- set r_

    for(label i = 0; i < nBins_; i++)
    {
       r_[i] = (0.5 + scalar(i)) * binWidth_;
    }

    //- set g_

    if (propsDict_.found("readCoefficients"))
    {
        readCoefficients_ = Switch(propsDict_.lookup("readCoefficients"));

        if(readCoefficients_)
        {
            dictionary coefficientsDict(propsDict_.subDict("coefficients"));
    
            setCoefficients(coefficientsDict);
        }
    }
    else
    {
        setCoefficients();
    }

    if (propsDict_.found("reducedUnits"))
    {
        reducedUnits_ = Switch(propsDict_.lookup("reducedUnits"));
    }

    createRdf();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

matteoliAndMansoori::~matteoliAndMansoori()
{}

void matteoliAndMansoori::setCoefficients(const dictionary& dict)
{
    h_ = readScalar(dict.lookup("h"));
    m_ = readScalar(dict.lookup("m"));
    gd_ = readScalar(dict.lookup("gd"));
    lambda_= readScalar(dict.lookup("lambda"));
    alpha_ = readScalar(dict.lookup("alpha"));
    beta_ = readScalar(dict.lookup("beta"));
    theta_ = readScalar(dict.lookup("theta"));

}

void matteoliAndMansoori::setCoefficients()
{
    h_ = (p1(403.5, -371.7, -1.552)/1000) + 1;
    gd_ = p1(1.708, -0.8569, 0.8196);
    lambda_ = p1(0.5644, -0.3057, 0.8579);
    beta_ = p1(5.289, -1.18, 0.3996);
    theta_ = p1(71.44, -46.68, 1.1);

    m_ = p2(22.79, -17.54, -0.0508);
    alpha_ = p2(0.2411, 0.1387, 4.216);
}


void matteoliAndMansoori::createRdf()
{
    for(label i = 0; i < g_.size(); i++)
    {
        scalar y = r_[i]/h_;

        if(y >= 1)
        {
            g_[i] = tailPart(y,m_,gd_,lambda_,alpha_,beta_);
        }
        else
        {
            g_[i] = initialPart(y,gd_,theta_);
        }
    }

    if(!reducedUnits_)
    {
        forAll(r_, r)
        {
            r_[r] *= 0.34e-9;
        }
    }
}

scalar matteoliAndMansoori::tailPart
(
    const scalar& y,
    const scalar& m,
    const scalar& gd,
    const scalar& lambda,
    const scalar& alpha,
    const scalar& beta
)
{
    return 1.0 + pow(y, -m)*(gd-1.0-lambda)
                    + (y - 1.0+lambda)*exp(-alpha*(y-1.0))*cos(beta*(y-1.0))/y;
}


scalar matteoliAndMansoori::initialPart
(

    const scalar& y,
    const scalar& gd,
    const scalar& theta
)
{
    return gd*exp(-theta*sqr(y-1));
}

scalar matteoliAndMansoori::p1
(
    const scalar& q1,
    const scalar& q2,
    const scalar& q3
)
{
    return (q1 + q2*exp(-1/T_))*exp(q3*density_);
}

scalar matteoliAndMansoori::p2
(
    const scalar& q1,
    const scalar& q2,
    const scalar& q3
)
{
    return (q1 + q2*exp(-1/T_))*((density_*density_)+q3)/density_;
}



// const scalarField& matteoliAndMansoori::g() const
// {
//     return g_;
// }
// 
// 
// const scalarField& matteoliAndMansoori::r() const
// {
//     return r_;
// }
// 
// const scalar& matteoliAndMansoori::binWidth() const
// {
//     return binWidth_;
// }



void matteoliAndMansoori::setRDF(radialDistribution& rdf, const Time& runTime)
{
    rdf.setRdf
    (
        name_,
        g_,
        r_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





} // End namespace Foam

// ************************************************************************* //
