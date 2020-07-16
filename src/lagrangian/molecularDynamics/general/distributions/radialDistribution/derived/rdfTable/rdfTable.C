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

Description

\*---------------------------------------------------------------------------*/

#include "rdfTable.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(rdfTable, 0);

addToRunTimeSelectionTable(rdfModel, rdfTable, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
rdfTable::rdfTable
(
//     Time& t,
//     const polyMesh& mesh,
//     moleculeCloud& molCloud,
    const dictionary& dict
)
:
    rdfModel(/*t, mesh, molCloud,*/ dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    name_(propsDict_.lookup("distributionName"))
//     rdf_()
{
//     word rdfName(propsDict_.lookup("name"));
//
//     readFromFile_ = Switch(propsDict_.lookup("readFromFile"));
//
//     if(!readFromFile_)
//     {
//         autoPtr<rdfModel> rdfTemp
//         (
//             rdfModel::New(t, propsDict_)
//         );
//
//         const scalarField& g = rdfTemp->g();
//         const scalarField& r = rdfTemp->r();
//         const scalar& binWidth = rdfTemp->binWidth();
//
//         rdf_ = radialDistribution(g, r, binWidth, rdfName);
//     }
//     else
//     {
//         rdf_ = radialDistribution(rdfName);
//
//         rdf_.readRDF(time_.time());
//     }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rdfTable::~rdfTable()
{}



void rdfTable::setRDF(radialDistribution& rdf, const Time& runTime)
{
    rdf.setRdfName(name_);
    rdf.readRDF(runTime);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





} // End namespace Foam

// ************************************************************************* //
