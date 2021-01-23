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

#include "Specie2Mixture.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::Specie2Mixture<MixtureType>::Specie2Mixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    MixtureType
    (
        thermoDict,
        mesh
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::nMoles
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).nMoles();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::W
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).W();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::W
(
    const word& specieName
) const
{
    return this->getLocalThermo(this->species()[specieName]).W();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::R
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).R();
}


template<class MixtureType>
Foam::label Foam::Specie2Mixture<MixtureType>::particleType
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).particleType();
}


template<class MixtureType>
Foam::label Foam::Specie2Mixture<MixtureType>::particleCharge
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).particleCharge();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::diameter
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).diameter();
}


template<class MixtureType>
Foam::DynamicList<Foam::scalar>
Foam::Specie2Mixture<MixtureType>::vibrationalList
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).vibrationalList();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::dissociationPotential
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).dissociationPotential();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::omega
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).omega();
}


template<class MixtureType>
Foam::label Foam::Specie2Mixture<MixtureType>::noVibrationalTemp
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).noVibrationalTemp();
}


template<class MixtureType>
Foam::label Foam::Specie2Mixture<MixtureType>::noElectronicLevels
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).noElectronicLevels();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::iHat
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).iHat();
}

// * * * * * * * * * * * * * Backward compatibility  * * * * * * * * * * * * //
template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cp
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Cp(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cv
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Cv(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Ha
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Ha(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hs
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Hs(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hc
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).Hc();
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::S
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).S(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Es
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Es(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::G
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).G(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::A
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).A(p, T);
}
// * * * * * * * * * * * * End backward compatibility  * * * * * * * * * * * //


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cp
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cp(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cv
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cv(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cv_t
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Cv_t(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cv_v
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cv_v(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cv_el
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cv_el(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cv_vel
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cv_vel(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cp_t
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Cp_t(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cp_v
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cp_v(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cp_el
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Cp_el(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Cp_vel
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    const scalar Cvvel = this->getLocalThermo(speciei).Cv_vel(p, Tve);
    
    scalar Cpvel = Cvvel;
    
    if (this->isElectron(speciei))
    {
        Cpvel += R(speciei);
    }
    
    return Cpvel;
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Ha
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Ha(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hs
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Hs(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hts
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Hts(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hvs
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Hvs(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hels
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Hels(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Hvels
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    const scalar Evel = this->getLocalThermo(speciei).Evels(p, Tve);
    
    scalar Hvel = Evel;
    
    if (this->isElectron(speciei))
    {
        Hvel += R(speciei)*Tve;
    }
    
    return Hvel;
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::S
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).S(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Es
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Es(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Ets
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Ets(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Evs
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Evs(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Eels
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Eels(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::Evels
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).Evels(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::HEt
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).HEt(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::HEv
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).HEv(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::HEel
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).HEel(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::HEvel
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).HEvel(p, Tve);
}


/*template<class MixtureType> TODO ABORTIVE WORK
Foam::scalar Foam::Specie2Mixture<MixtureType>::HEvel_mode
(
    const label speciei,
    const label mode,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).HEvel_mode(mode, p, Tve);
}*/


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::Specie2Mixture<MixtureType>::Cv_t
(
    const label speciei,
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCvt(new scalarField(T.size()));
    scalarField& Cvt = tCvt.ref();

    forAll(T, facei)
    {
        Cvt[facei] = Cv_t(speciei, p[facei], T[facei]);
    }

    return tCvt;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::Specie2Mixture<MixtureType>::Cv_vel
(
    const label speciei,
    const scalarField& p,
    const scalarField& Tve,
    const label patchi
) const
{
    tmp<scalarField> tCvvel(new scalarField(Tve.size()));
    scalarField& Cvvel = tCvvel.ref();

    forAll(Tve, facei)
    {
        Cvvel[facei] = Cv_vel(speciei, p[facei], Tve[facei]);
    }

    return tCvvel;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::Specie2Mixture<MixtureType>::het
(
    const label speciei,
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> thet(new scalarField(T.size()));
    scalarField& het = thet.ref();

    forAll(T, celli)
    {
        het[celli] = HEt(speciei, p[celli], T[celli]);
    }

    return thet;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::Specie2Mixture<MixtureType>::het
(
    const label speciei,
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> thet(new scalarField(T.size()));
    scalarField& het = thet.ref();

    forAll(T, facei)
    {
        het[facei] = HEt(speciei, p[facei], T[facei]);
    }

    return thet;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::Specie2Mixture<MixtureType>::hevel
(
    const label speciei,
    const scalarField& p,
    const scalarField& Tve,
    const labelList& cells
) const
{
    tmp<scalarField> thevel(new scalarField(Tve.size()));
    scalarField& hevel = thevel.ref();

    forAll(Tve, celli)
    {
        hevel[celli] = HEvel(speciei, p[celli], Tve[celli]);
    }

    return thevel;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::Specie2Mixture<MixtureType>::hevel
(
    const label speciei,
    const scalarField& p,
    const scalarField& Tve,
    const label patchi
) const
{
    tmp<scalarField> thevel(new scalarField(Tve.size()));
    scalarField& hevel = thevel.ref();

    forAll(Tve, facei)
    {
        hevel[facei] = HEvel(speciei, p[facei], Tve[facei]);
    }

    return thevel;
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::TtHEt
(
    const label speciei,
    const scalar het,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).TtHEt(het, p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::zetar
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).zetar(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::zetav
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).zetav(p, Tve);
}


/*template<class MixtureType> TODO ABORTIVE WORK
Foam::scalar Foam::Specie2Mixture<MixtureType>::zetav_mode
(
    const label speciei,
    const label mode,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).zetav_mode(mode, p, Tve);
}*/


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::zetael
(
    const label speciei,
    const scalar p,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).zetael(p, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::G
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).G(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::A
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).A(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::psi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).psi(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::mu
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).mu(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::kappatr
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).kappatr(p, T);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::kappave
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return this->getLocalThermo(speciei).kappave(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::kappa
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    return
        this->getLocalThermo(speciei).kappatr(p, T)
      + this->getLocalThermo(speciei).kappave(p, T, Tve);
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::alphatr
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    const scalar Cptr = max(this->getLocalThermo(speciei).Cp_t(p, T), 1e-12);
    
    return this->getLocalThermo(speciei).kappatr(p, T)/Cptr;
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::alphave
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    const scalar loc_Cpvel = max(Cp_vel(speciei, p, Tve), 1e-12);
    
    return this->getLocalThermo(speciei).kappave(p, T, Tve)/loc_Cpvel;
}


template<class MixtureType>
Foam::scalar Foam::Specie2Mixture<MixtureType>::alpha
(
    const label speciei,
    const scalar p,
    const scalar T,
    const scalar Tve
) const
{
    const scalar loc_Cp = max(Cp(speciei, p, T, Tve), 1e-12);
    
    return kappa(speciei, p, T, Tve)/loc_Cp;
}


// ************************************************************************* //
