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

\*---------------------------------------------------------------------------*/

#include "mhdModel.H"
#include "electricalConductivityModel.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(mhdModel, 0);
        defineRunTimeSelectionTable(mhdModel, thermo);
        defineRunTimeSelectionTable(mhdModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

IOobject mhdModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "mhdProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOobject>())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void mhdModel::initialise()
{
    if (active_)
    {
        Info << "Initialising MHD model" << endl;
        
        electricalConductivity_.reset
        (
            electricalConductivityModel::New(*this, mesh_).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mhdModel::mhdModel(const rho2ReactionThermo& thermo)
:
    IOdictionary
    (
        IOobject
        (
            "mhdProperties",
            thermo.T().time().constant(),
            thermo.T().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(thermo.T().mesh()),
    time_(thermo.T().time()),
    thermo_(thermo),
    active_(true),
    hallEffect_(false),
    coeffs_(dictionary::null),
    electricalConductivity_(NULL)
{
    const dictionary mhdProperties =
    (
        IFstream
        (
            fileName(thermo.lookup("mhdDictFile")).expand()
        )()
    );
        
    *this <<= mhdProperties;
    
    initialise();
    
    hallEffect_ = lookupOrDefault("hallEffect", false);
}


mhdModel::mhdModel
(
    const word& type,
    const rho2ReactionThermo& thermo
)
:
    IOdictionary(createIOobject(thermo.T().mesh())),
    mesh_(thermo.T().mesh()),
    time_(thermo.T().time()),
    thermo_(thermo),
    active_(lookupOrDefault("active", true)),
    hallEffect_(lookupOrDefault("hallEffect", false)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    electricalConductivity_(NULL)
{
    if (readOpt() == IOobject::NO_READ)
    {
        active_ = false;
    }
    
    initialise();
}


mhdModel::mhdModel
(
    const word& type,
    const dictionary& dict,
    const rho2ReactionThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            "mhdProperties",
            thermo.T().time().constant(),
            thermo.T().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dict
    ),
    mesh_(thermo.T().mesh()),
    time_(thermo.T().time()),
    thermo_(thermo),
    active_(lookupOrDefault("active", true)),
    hallEffect_(lookupOrDefault("hallEffect", false)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    electricalConductivity_(NULL)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

mhdModel::~mhdModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool mhdModel::read()
{
    if (regIOobject::read())
    {
        lookup("active") >> active_;
        coeffs_ = subOrEmptyDict(type() + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}


const electricalConductivityModel&
mhdModel::electricalConductivity() const
{
    if (!electricalConductivity_.valid())
    {
        FatalErrorIn
        (
            "const Foam::mhd::electricalConductivityModel&"
            "Foam::mhd::mhdModel::electricalConductivity() const"
        )   << "Requested electrical conductivity model, but model is "
            << "not activated" << abort(FatalError);
    }

    return electricalConductivity_();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam


// ************************************************************************* //
