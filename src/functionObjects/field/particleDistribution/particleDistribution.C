/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "particleDistribution.H"
#include "addToRunTimeSelectionTable.H"
#include "general.H"
#include "fvMesh.H"
#include "cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(particleDistribution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        particleDistribution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::particleDistribution::particleDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(runTime, name),
    cloudName_("unknown-cloudName"),
    nameVsBinWidth_(),
    tagFieldName_("none"),
    rndGen_(1234, -1),
    writerPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::particleDistribution::~particleDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::particleDistribution::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.lookup("cloud") >> cloudName_;
        dict.lookup("nameVsBinWidth") >> nameVsBinWidth_;
        dict.readIfPresent("tagField", tagFieldName_);
        word format(dict.lookup("setFormat"));
        writerPtr_ = writer<scalar>::New(format);

        Info<< type() << " " << name() << " output:" << nl
            << "    Processing cloud : " << cloudName_ << nl
            << endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::particleDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::particleDistribution::write()
{
    Log << type() << " " << name() << " output:" << endl;

    if (!mesh_.foundObject<cloud>(cloudName_))
    {
        wordList cloudNames(mesh_.names<cloud>());

        WarningInFunction
            << "Unable to find cloud " << cloudName_
            << " in the mesh database.  Available clouds include:"
            << cloudNames << endl;

        return false;
    }

    const cloud& c = mesh_.lookupObject<cloud>(cloudName_);

    objectRegistry cloudObr
    (
        IOobject
        (
            name() & "CloudRegistry",
            mesh_.time().timeName(),
            cloud::prefix,
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    c.writeObjects(cloudObr);

    List<DynamicList<label>> tagAddr;
    if
    (
        tagFieldName_ != "none"
     && cloudObr.foundObject<IOField<scalar>>(tagFieldName_)
    )
    {
        // Tag field present - generate distribution per tag
        const IOField<label>& tag =
            cloudObr.lookupObject<IOField<label>>(tagFieldName_);
        const HashSet<label> tagMap(tag);
        const label tagMax = tagMap.size();

        List<DynamicList<label>> tagAddr(tagMax);
        forAll(tag, i)
        {
            label newTag = tagMap[tag[i]];
            tagAddr[newTag].append(i);
        }
    }


    bool ok = false;
    forAll(nameVsBinWidth_, i)
    {
        ok = false;
        ok = ok || processField<scalar>(cloudObr, i, tagAddr);
        ok = ok || processField<vector>(cloudObr, i, tagAddr);
        ok = ok || processField<tensor>(cloudObr, i, tagAddr);
        ok = ok || processField<sphericalTensor>(cloudObr, i, tagAddr);
        ok = ok || processField<symmTensor>(cloudObr, i, tagAddr);
        ok = ok || processField<tensor>(cloudObr, i, tagAddr);

        if (log && !ok)
        {
            WarningInFunction
                << "Unable to find field " << nameVsBinWidth_[i].first()
                << " in the " << cloudName_ << " cloud database" << endl;
        }
    }

    return true;
}


void Foam::functionObjects::particleDistribution::generateDistribution
(
    const word& fieldName,
    const scalarField& field,
    const scalar binWidth,
    const label tag
)
{
    if (field.empty())
    {
        return;
    }

    word fName(fieldName);
    if (tag != -1)
    {
        fName = fName + '_' + Foam::name(tag);
    }

    distributionModels::general distribution
    (
        field,
        binWidth,
        rndGen_
    );

    const Field<scalar> distX(distribution.x());
    const Field<scalar> distY(distribution.y());

    pointField xBin(distX.size(), Zero);
    xBin.replace(0, distX);
    const coordSet coords
    (
        fName,
        "x",
        xBin,
        distX
    );

    const wordList fieldNames(1, fName);

    fileName outputPath(baseTimeDir());
    mkDir(outputPath);
    OFstream graphFile(outputPath/writerPtr_->getFileName(coords, fieldNames));

    Log << "    Writing distribution of " << fieldName
        << " to " << graphFile.name() << endl;

    List<const scalarField*> yPtrs(1);
    yPtrs[0] = &distY;
    writerPtr_->write(coords, fieldNames, yPtrs, graphFile);
}


// ************************************************************************* //
