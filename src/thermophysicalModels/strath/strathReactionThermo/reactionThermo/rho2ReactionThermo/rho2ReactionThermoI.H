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

#include "rho2ReactionThermo.H"
#include "PtrList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::rho2ReactionThermo::Temperature_NewtonRaphson
(
    scalar f,
    scalar p,
    scalar T0,
    const PtrList<scalar>& YList,
    scalar (Foam::basic2MultiComponentMixture::*F)
    (
        const label,
        const scalar,
        const scalar
    ) const,
    scalar (Foam::basic2MultiComponentMixture::*dFdT)
    (
        const label,
        const scalar,
        const scalar
    ) const,
    scalar (Foam::rho2ReactionThermo::*limit)(const scalar) const
)
{
    // The tolerance and maximum number of iterations are identical to the ones
    // in multiThermo.C
    scalar tol_ = 1e-6;
    label maxIter_ = 100;

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    label iter = 0;

    do
    {
        Test = Tnew;

        scalar Fipt(0), dFipt(0);
        forAll(composition().species(), i)
        {
            Fipt  += (composition().*F)(i, p, Test)*YList[i];
            dFipt += (composition().*dFdT)(i, p, Test)*YList[i];
        }

        Tnew =
            (this->*limit)
            (Test - (Fipt - f)/dFipt);

        if (iter++ > maxIter_)
        {
            // the temperature cannot be derived from the energy
            // it remains unchanged
            Tnew = T0;
            Test = T0;
            Foam::rho2ReactionThermo::hasCrashedButRecovered = true;
        }

    } while (mag(Tnew - Test) > Ttol);
    
    return Tnew;
}


Foam::scalar Foam::rho2ReactionThermo::Temperature_Secant
(
    scalar f,
    scalar p,
    scalar T0,
    const scalar Tt,
    const PtrList<scalar>& YList,
    scalar (Foam::basic2MultiComponentMixture::*F)
    (
        const label,
        const scalar,
        const scalar
    ) const,
    scalar (Foam::rho2ReactionThermo::*limit)(const scalar) const
)
{
    // The tolerance and maximum number of iterations are identical to the ones
    // in multiThermo.C
    scalar tol_ = 1e-6;
    label maxIter_ = 100;

    scalar Test2 = Tt;
    scalar Test1 = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    label iter = 0;
    
    do
    {
        Test1 = Tnew;

        scalar Fipt1(0), Fipt2(0);
        forAll(composition().species(), i)
        {
            Fipt1 += (composition().*F)(i, p, Test1)*YList[i];
            Fipt2 += (composition().*F)(i, p, Test2)*YList[i];
        }

        Tnew =
            (this->*limit)
            (Test1 - (Fipt1 - f)*(Test1 - Test2)/(Fipt1 - Fipt2));

        if (iter++ > maxIter_)
        {
            // the temperature cannot be derived from the energy
            // it remains unchanged
            Tnew = T0;
            Test1 = T0;
            Foam::rho2ReactionThermo::hasCrashedButRecovered = true;
        }
        
/*        Info << iter << tab << Test2 << tab << Test1 << tab << Tnew << endl; */
        
        Test2 = Test1;

    } while (mag(Tnew - Test1) > Ttol);

    return Tnew;
}


Foam::scalar Foam::rho2ReactionThermo::Temperature
(
    label i,
    scalar f,
    scalar p,
    scalar T0,
    scalar (Foam::basic2MultiComponentMixture::*F)
    (
        const label,
        const scalar,
        const scalar
    ) const,
    scalar (Foam::basic2MultiComponentMixture::*dFdT)
    (
        const label,
        const scalar,
        const scalar
    ) const,
    scalar (Foam::rho2ReactionThermo::*limit)(const scalar) const
)
{
    // The tolerance and maximum number of iterations are identical to the ones
    // in multiThermo.C
    scalar tol_ = 1e-6;
    label maxIter_ = 100;

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    label iter = 0;

    do
    {
        Test = Tnew;

        scalar Fipt((composition().*F)(i, p, Test)),
            dFipt((composition().*dFdT)(i, p, Test));

        Tnew =
            (this->*limit)
            (Test - (Fipt - f)/dFipt);

        if (iter++ > maxIter_)
        {
            // the temperature cannot be derived from the energy
            // it remains unchanged
            Tnew = T0;
            Test = T0;
            Foam::rho2ReactionThermo::hasCrashedButRecovered = true;
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


Foam::scalar Foam::rho2ReactionThermo::TEs
(
    const scalar e,
    const scalar p,
    const scalar T0,
    const PtrList<scalar>& YList
)
{
    return Temperature_NewtonRaphson
    (
        e,
        p,
        T0,
        YList,
        &Foam::basic2MultiComponentMixture::Es,
        &Foam::basic2MultiComponentMixture::Cv,
        &Foam::rho2ReactionThermo::limit
    );
}


Foam::scalar Foam::rho2ReactionThermo::TtEts
(
    const scalar e,
    const scalar p,
    const scalar T0,
    const PtrList<scalar>& YList
)
{
    return Temperature_NewtonRaphson
    (
        e,
        p,
        T0,
        YList,
        &Foam::basic2MultiComponentMixture::Ets,
        &Foam::basic2MultiComponentMixture::Cv_t,
        &Foam::rho2ReactionThermo::limit
    );
}


Foam::scalar Foam::rho2ReactionThermo::TvelEvels
(
    const scalar e,
    const scalar p,
    const scalar T0,
    const PtrList<scalar>& YList
)
{
    return Temperature_NewtonRaphson
    (
        e,
        p,
        T0,
        YList,
        &Foam::basic2MultiComponentMixture::Evels,
        &Foam::basic2MultiComponentMixture::Cv_vel,
        &Foam::rho2ReactionThermo::limit
    );
}


Foam::scalar Foam::rho2ReactionThermo::TvelEvels
(
    const scalar e,
    const scalar p,
    const scalar T0,
    const scalar Tt,
    const PtrList<scalar>& YList
)
{
    return Temperature_Secant
    (
        e,
        p,
        T0,
        Tt,
        YList,
        &Foam::basic2MultiComponentMixture::Evels,
        &Foam::rho2ReactionThermo::limit
    );
}


Foam::scalar Foam::rho2ReactionThermo::TvelEvels
(
    const label i,
    const scalar e,
    const scalar p,
    const scalar T0
)
{
    return Temperature
    (
        i,
        e,
        p,
        T0,
        &Foam::basic2MultiComponentMixture::Evels,
        &Foam::basic2MultiComponentMixture::Cv_vel,
        &Foam::rho2ReactionThermo::limit
    );
}


Foam::scalar Foam::rho2ReactionThermo::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        temperatureFieldOutOfRange = true;
        noCellsWithTemperatureFieldOutOfRange += 1;
        minMaxTemperatureFieldOutOfRange[0] = min(T, Tlow_);
        minMaxTemperatureFieldOutOfRange[1] = max(T, Thigh_);

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}

// ************************************************************************* //
