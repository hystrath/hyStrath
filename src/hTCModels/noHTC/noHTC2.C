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

#include "noHTC2.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class HTempThermoType>
Foam::hTC2Models::noHTC2<HTempThermoType>::noHTC2
(
    const word& modelType,
    const fvMesh& mesh
)
:
    HTempThermoType(modelType, mesh)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class HTempThermoType>
Foam::hTC2Models::noHTC2<HTempThermoType>::~noHTC2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class HTempThermoType>
void Foam::hTC2Models::noHTC2<HTempThermoType>::correct()
{
//  Do Nothing
}


template<class HTempThermoType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::hTC2Models::noHTC2<HTempThermoType>::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
}


template<class HTempThermoType>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::noHTC2<HTempThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
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

    return tdQ;
}


template<class HTempThermoType>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::noHTC2<HTempThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
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

    return tSh;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::noHTC2<Type>::Scv() const
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

    return tScv;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::noHTC2<Type>::Scv(const label i) const
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

    return tScv;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::noHTC2<Type>::Seiir() const
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

    return tSeiir;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::noHTC2<Type>::Seiir(const label i) const
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

    return tSeiir;
}


template<class HTempThermoType>
bool Foam::hTC2Models::noHTC2<HTempThermoType>::read()
{
    if (HTempThermoType::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
