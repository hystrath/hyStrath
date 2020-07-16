/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
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
    rdf

Description

\*----------------------------------------------------------------------------*/

#include "rdf.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// constructor
rdf::rdf
(
    const dictionary& dict,
    Time& t
)
:
    rdfDistr_(),
    rdfModel_(),
    time_(t)
{

    rdfModel_ = autoPtr<rdfModel>
    (
        rdfModel::New(dict)
    );

    rdfModel_->setRDF(rdfDistr_, time_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rdf::~rdf()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const radialDistribution& rdf::RDF() const
{
    return rdfDistr_;
}

void rdf::write()
{
    rdfDistr_.writeRDF(time_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
