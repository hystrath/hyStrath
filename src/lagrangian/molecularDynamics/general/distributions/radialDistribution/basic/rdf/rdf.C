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
