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

#include "laminar2.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::hTC2Models::laminar2<Type>::laminar2
(
    const word& modelType,
    const fvMesh& mesh
)
:
    Type(modelType, mesh),
    integrateReactionRate_
    (
        this->coeffs().lookupOrDefault("integrateReactionRate", true)
    )
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::hTC2Models::laminar2<Type>::~laminar2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class Type>
void Foam::hTC2Models::laminar2<Type>::correct()
{
    if (this->active())
    {
        if (integrateReactionRate_)
        {
            word ddtScheme(this->mesh().ddtScheme("Yi"));

            if (ddtScheme == fv::localEulerDdtScheme<scalar>::typeName)
            {
                const scalarField& rDeltaT =
                    this->mesh().objectRegistry::
                    template lookupObject<volScalarField>
                    (
                        "rDeltaT"
                    );

                if (this->coeffs().found("maxIntegrationTime"))
                {
                    scalar maxIntegrationTime
                    (
                        readScalar(this->coeffs().lookup("maxIntegrationTime"))
                    );

                    this->chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxIntegrationTime)()
                    );
                }
                else
                {
                    this->chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                this->chemistryPtr_->solve(this->mesh().time().deltaTValue());
            }
        }
        else
        {
            this->chemistryPtr_->calculate();
        }
    }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::hTC2Models::laminar2<Type>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    if (this->active())
    {
        const label specieI = this->thermo().composition().species()[Y.name()];

        Su += this->chemistryPtr_->RR(specieI);
    }
    return tSu;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
        )
    );

    if (this->active())
    {
        tdQ.ref() = this->chemistryPtr_->dQ();
    }

    return tdQ;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    if (this->active())
    {
        tSh.ref() = this->chemistryPtr_->Sh();
    }

    return tSh;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::Scv() const
{
    tmp<volScalarField> tScv
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Scv",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    if (this->active())
    {
        tScv.ref() = this->chemistryPtr_->Scv();
    }

    return tScv;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::Scv(const label i) const
{
    tmp<volScalarField> tScv
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Scv_" + word(i),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    if (this->active())
    {
        tScv.ref() = this->chemistryPtr_->Scv(i);
    }

    return tScv;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::Seiir() const
{
    tmp<volScalarField> tSeiir
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Seiir",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    if (this->active())
    {
        tSeiir.ref() = this->chemistryPtr_->Seiir();
    }

    return tSeiir;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::laminar2<Type>::Seiir(const label i) const
{
    tmp<volScalarField> tSeiir
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Seiir_" + word(i),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    if (this->active())
    {
        tSeiir.ref() = this->chemistryPtr_->Seiir(i);
    }

    return tSeiir;
}


template<class Type>
bool Foam::hTC2Models::laminar2<Type>::read()
{
    if (Type::read())
    {
        /*this->coeffs().lookupOrDefault("integrateReactionRate", true)
            >> integrateReactionRate_;*/ // DELETED VINCENT 07/09/2016
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
