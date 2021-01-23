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

#include "rho2Thermo.H"

#include "zeroGradientFvPatchFields.H" // NEW VINCENT 13/04/2016
#include "fixedRhoFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rho2Thermo, 0);
    defineRunTimeSelectionTable(rho2Thermo, fvMesh);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::rho2Thermo::rhoBoundaryTypes() // NEW VINCENT 13/04/2016
{
    const volScalarField::Boundary& pbf = this->T_.boundaryField();
    wordList rhoBoundaryTypes = pbf.types();

    forAll(rhoBoundaryTypes, patchi)
    {
        if (rhoBoundaryTypes[patchi] == "waveTransmissive")
        {
            rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
        }
        else if (pbf[patchi].fixesValue())
        {
            rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField::typeName;
        }
    }

    return rhoBoundaryTypes;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rho2Thermo::rho2Thermo(const fvMesh& mesh, const word& phaseName)
:
    multi2Thermo(mesh, phaseName),

    rho_
    (
        IOobject
        (
            phasePropertyName("rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity,
        rhoBoundaryTypes()
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 1e-6)
    ),
    
    kappatr_
    (
        IOobject
        (
            phasePropertyName("kappatr"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0)
    ),

    kappave_
    (
        IOobject
        (
            phasePropertyName("kappave"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0)
    ),
    
    kappa_
    (
        IOobject
        (
            phasePropertyName("kappa"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0)
    ),

    alphatr_
    (
        IOobject
        (
            phasePropertyName("alphatr"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimLength/dimTime
    ),

    alphave_
    (
        IOobject
        (
            phasePropertyName("alphave"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimLength/dimTime
    ),
    
    alpha_
    (
        IOobject
        (
            phasePropertyName("alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimLength/dimTime
    )
{}


Foam::rho2Thermo::rho2Thermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    multi2Thermo(mesh, dict, phaseName),

    rho_
    (
        IOobject
        (
            phasePropertyName("rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity,
        rhoBoundaryTypes()
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 1e-6)
    ),

    kappatr_
    (
        IOobject
        (
            phasePropertyName("kappatr"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0)
    ),

    kappave_
    (
        IOobject
        (
            phasePropertyName("kappave"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0)
    ),
    
    kappa_
    (
        IOobject
        (
            phasePropertyName("kappa"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0)
    ),

    alphatr_
    (
        IOobject
        (
            phasePropertyName("alphatr"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimLength/dimTime
    ),

    alphave_
    (
        IOobject
        (
            phasePropertyName("alphave"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimLength/dimTime
    ),
    
    alpha_
    (
        IOobject
        (
            phasePropertyName("alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimLength/dimTime
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rho2Thermo> Foam::rho2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basic2Thermo::New<rho2Thermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rho2Thermo::~rho2Thermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rho2Thermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::rho2Thermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::rho()
{
    return rho_;
}


const Foam::volScalarField& Foam::rho2Thermo::psi() const
{
    return psi_;
}


Foam::volScalarField& Foam::rho2Thermo::psi()
{
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::rho2Thermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::rho2Thermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::mu()
{
    return mu_;
}


const Foam::volScalarField& Foam::rho2Thermo::kappatr() const
{
    return kappatr_;
}


const Foam::scalarField& Foam::rho2Thermo::kappatr(const label patchi) const
{
    return kappatr_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::kappatr()
{
    return kappatr_;
}


const Foam::volScalarField& Foam::rho2Thermo::kappave() const
{
    return kappave_;
}


const Foam::scalarField& Foam::rho2Thermo::kappave(const label patchi) const
{
    return kappave_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::kappave()
{
    return kappave_;
}


const Foam::volScalarField& Foam::rho2Thermo::kappa() const
{
    return kappa_;
}


const Foam::scalarField& Foam::rho2Thermo::kappa(const label patchi) const
{
    return kappa_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::kappa()
{
    return kappa_;
}


const Foam::volScalarField& Foam::rho2Thermo::alphatr() const
{
    return alphatr_;
}


const Foam::scalarField& Foam::rho2Thermo::alphatr(const label patchi) const
{
    return alphatr_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::alphatr()
{
    return alphatr_;
}


const Foam::volScalarField& Foam::rho2Thermo::alphave() const
{
    return alphave_;
}


const Foam::scalarField& Foam::rho2Thermo::alphave(const label patchi) const
{
    return alphave_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::alphave()
{
    return alphave_;
}


const Foam::volScalarField& Foam::rho2Thermo::alpha() const
{
    return alpha_;
}


const Foam::scalarField& Foam::rho2Thermo::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rho2Thermo::alpha()
{
    return alpha_;
}


Foam::tmp<Foam::volScalarField> Foam::rho2Thermo::kappaEff
(
    const volScalarField& alphat
) const
{
    const fvMesh& mesh = alphat.mesh();
    
    tmp<volScalarField> tkappaEff
    (
        new volScalarField
        (
            IOobject
            (
                "kappaEff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimensionSet(1, 1, -3, -1, 0),
                0.0
            )
        )
    );
    
    volScalarField& kappaEff = tkappaEff.ref();
    
    kappaEff = kappa() + CpMix()*alphat;
      
    return tkappaEff;
}


Foam::tmp<Foam::scalarField> Foam::rho2Thermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return kappa(patchi)+ CpMix(patchi)*alphat;
}


Foam::tmp<Foam::volScalarField> Foam::rho2Thermo::alphaEff
(
    const volScalarField& alphat
) const
{
    const fvMesh& mesh = alphat.mesh();
    
    tmp<volScalarField> talphaEff
    (
        new volScalarField
        (
            IOobject
            (
                "alphaEff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimMass/dimLength/dimTime,
                0.0
            )
        )
    );
    
    volScalarField& alphaEff = talphaEff.ref();
    
    alphaEff = alpha() + alphat;
    
    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::rho2Thermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return alphatr(patchi) + alphat;
}


// ************************************************************************* //
