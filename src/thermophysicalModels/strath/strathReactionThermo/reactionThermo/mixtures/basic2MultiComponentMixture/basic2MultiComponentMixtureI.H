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

inline const Foam::label&
Foam::basic2MultiComponentMixture::vibTempAssociativity
(
    const label speciei
) const // NEW VINCENT 05/03/2016
{
    return vibTempAssociativity_[speciei];
}


inline const Foam::label&
Foam::basic2MultiComponentMixture::vibTempAssociativity
(
    const word& specieName
) const // NEW VINCENT 06/08/2016
{
    return vibTempAssociativity_[species_[specieName]];
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::e()
{
    return *e_;
}


inline const Foam::volScalarField&
Foam::basic2MultiComponentMixture::e() const
{
    return *e_;
}


//----                               Y                                  -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}


//----                               X                                  -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::X()
{
    return X_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::X() const
{
    return X_;
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::X(const label i)
{
    return X_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::X
(
    const label i
) const
{
    return X_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::X
(
    const word& specieName
)
{
    return X_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::X
(
    const word& specieName
) const
{
    return X_[species_[specieName]];
}


//----                               nD                                 -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::nD()
{
    return nD_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::nD() const
{
    return nD_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::nD(const label i)
{
    return nD_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::nD
(
    const label i
) const
{
    return nD_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::nD
(
    const word& specieName
)
{
    return nD_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::nD
(
    const word& specieName
) const
{
    return nD_[species_[specieName]];
}


//----                               pP                                 -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::pP()
{
    return pP_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::pP() const
{
    return pP_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::pP(const label i)
{
    return pP_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::pP
(
    const label i
) const
{
    return pP_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::pP
(
    const word& specieName
)
{
    return pP_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::pP
(
    const word& specieName
) const
{
    return pP_[species_[specieName]];
}


//----                               pD                                 -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::pD()
{
    return pD_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::pD() const
{
    return pD_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::pD(const label i)
{
    return pD_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::pD
(
    const label i
) const
{
    return pD_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::pD
(
    const word& specieName
)
{
    return pD_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::pD
(
    const word& specieName
) const
{
    return pD_[species_[specieName]];
}


//----                               hev                                -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::hev()
{
    return hev_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::hev() const
{
    return hev_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::hev(const label i)
{
    return hev_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::hev
(
    const label i
) const
{
    return hev_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::hev
(
    const word& specieName
)
{
    return hev_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::hev
(
    const word& specieName
) const
{
    return hev_[species_[specieName]];
}


//----                               heel                               -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::heel()
{
    return heel_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::heel() const
{
    return heel_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::heel(const label i)
{
    return heel_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::heel
(
    const label i
) const
{
    return heel_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::heel
(
    const word& specieName
)
{
    return heel_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::heel
(
    const word& specieName
) const
{
    return heel_[species_[specieName]];
}


//----                               hevel                              -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::hevel()
{
    return hevel_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::hevel() const
{
    return hevel_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::hevel(const label i)
{
    return hevel_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel
(
    const label i
) const
{
    return hevel_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel
(
    const word& specieName
)
{
    return hevel_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel
(
    const word& specieName
) const
{
    return hevel_[species_[specieName]];
}


//----                               h                                  -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::enthalpy()
{
    return h_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::enthalpy() const
{
    return h_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::enthalpy(const label i)
{
    return h_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::enthalpy
(
    const label i
) const
{
    return h_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::enthalpy
(
    const word& specieName
)
{
    return h_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::enthalpy
(
    const word& specieName
) const
{
    return h_[species_[specieName]];
}


//----                               hevel_mode                         -----// TODO ABORTIVE WORK
/*inline Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::basic2MultiComponentMixture::hevel_mode()
{
    return modehevel_;
}


inline const Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::basic2MultiComponentMixture::hevel_mode() const
{
    return modehevel_;
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel_mode(const label i, const label m)
{
    return modehevel_[i][m];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel_mode
(
    const label i, const label m
) const
{
    return modehevel_[i][m];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel_mode
(
    const word& specieName, const label m
)
{
    return modehevel_[species_[specieName]][m];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::hevel_mode
(
    const word& specieName, const label m
) const
{
    return modehevel_[species_[specieName]][m];
}*/


//----                               Tv                                 -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::Tv()
{
    return spTv_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::Tv() const
{
    return spTv_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::Tv(const label i)
{
    return spTv_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv
(
    const label i
) const
{
    return spTv_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv
(
    const word& specieName
)
{
    return spTv_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv
(
    const word& specieName
) const
{
    return spTv_[species_[specieName]];
}


//----                               Tv_mode                            -----// TODO ABORTIVE WORK
/*inline Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::basic2MultiComponentMixture::Tv_mode()
{
    return spmodeTv_;
}


inline const Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::basic2MultiComponentMixture::Tv_mode() const
{
    return spmodeTv_;
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv_mode(const label i, const label m)
{
    return spmodeTv_[i][m];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv_mode
(
    const label i, const label m
) const
{
    return spmodeTv_[i][m];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv_mode
(
    const word& specieName, const label m
)
{
    return spmodeTv_[species_[specieName]][m];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::Tv_mode
(
    const word& specieName, const label m
) const
{
    return spmodeTv_[species_[specieName]][m];
}*/


//----                               zetar                              -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::zetar()
{
    return zetaRot_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::zetar() const
{
    return zetaRot_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::zetar(const label i)
{
    return zetaRot_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetar
(
    const label i
) const
{
    return zetaRot_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::zetar
(
    const word& specieName
)
{
    return zetaRot_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetar
(
    const word& specieName
) const
{
    return zetaRot_[species_[specieName]];
}


//----                               zetav                              -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::zetav()
{
    return zetaVib_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::zetav() const
{
    return zetaVib_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::zetav(const label i)
{
    return zetaVib_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav
(
    const label i
) const
{
    return zetaVib_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav
(
    const word& specieName
)
{
    return zetaVib_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav
(
    const word& specieName
) const
{
    return zetaVib_[species_[specieName]];
}


//----                               zetav_mode                         -----// TODO ABORTIVE WORK
/*inline Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::basic2MultiComponentMixture::zetav_mode()
{
    return modezetaVib_;
}


inline const Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::basic2MultiComponentMixture::zetav_mode() const
{
    return modezetaVib_;
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav_mode(const label i, const label m)
{
    return modezetaVib_[i][m];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav_mode
(
    const label i, const label m
) const
{
    return modezetaVib_[i][m];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav_mode
(
    const word& specieName, const label m
)
{
    return modezetaVib_[species_[specieName]][m];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetav_mode
(
    const word& specieName, const label m
) const
{
    return modezetaVib_[species_[specieName]][m];
}*/


//----                               zetael                             -----//
inline Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::zetael()
{
    return zetaElec_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basic2MultiComponentMixture::zetael() const
{
    return zetaElec_;
}


inline Foam::volScalarField&
Foam::basic2MultiComponentMixture::zetael(const label i)
{
    return zetaElec_[i];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetael
(
    const label i
) const
{
    return zetaElec_[i];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::zetael
(
    const word& specieName
)
{
    return zetaElec_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basic2MultiComponentMixture::zetael
(
    const word& specieName
) const
{
    return zetaElec_[species_[specieName]];
}


inline Foam::volScalarField& Foam::basic2MultiComponentMixture::Wmix()
{
    return Wmix_;
}


inline const Foam::volScalarField&
Foam::basic2MultiComponentMixture::Wmix() const
{
    return Wmix_;
}


inline bool Foam::basic2MultiComponentMixture::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}


inline bool Foam::basic2MultiComponentMixture::isMolecule
(
    const word& specieName
) const
{
    return molecules_.contains(specieName);
}


inline bool Foam::basic2MultiComponentMixture::isMolecule
(
    const label speciei
) const
{
    forAll(moleculeIds_, i)
    {
        if (speciei == moleculeIds_[i])
        {
            return true;
        }
    }
    
    return false;
}


inline bool Foam::basic2MultiComponentMixture::isAtom
(
    const word& specieName
) const
{
    return atoms_.contains(specieName);
}


inline bool Foam::basic2MultiComponentMixture::isAtom
(
    const label speciei
) const
{
    forAll(atomIds_, i)
    {
        if (speciei == atomIds_[i])
        {
            return true;
        }
    }
    
    return false;
}


inline bool Foam::basic2MultiComponentMixture::isIon
(
    const word& specieName
) const
{
    return ions_.contains(specieName);
}


inline bool Foam::basic2MultiComponentMixture::isIon
(
    const label speciei
) const
{
    forAll(ionIds_, i)
    {
        if (speciei == ionIds_[i])
        {
            return true;
        }
    }
    
    return false;
}


inline bool Foam::basic2MultiComponentMixture::containsVib
(
    const word& specieName
) const
{
    return solvedVibEqSpecies_.contains(specieName);
}


inline bool Foam::basic2MultiComponentMixture::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>& Foam::basic2MultiComponentMixture::active() const
{
    return active_;
}


inline void Foam::basic2MultiComponentMixture::setActive(label speciei)
{
    active_[speciei] = true;
}


inline void Foam::basic2MultiComponentMixture::setInactive(label speciei)
{
    active_[speciei] = false;
}


inline Foam::scalar Foam::basic2MultiComponentMixture::fres
(
    const scalar ft,
    const scalar stoicRatio
) const
{
    return max(ft - (1.0 - ft)/stoicRatio, 0.0);
}


inline Foam::tmp<Foam::volScalarField> Foam::basic2MultiComponentMixture::fres
(
    const volScalarField& ft,
    const dimensionedScalar& stoicRatio
) const
{
    return max(ft - (scalar(1) - ft)/stoicRatio.value(), scalar(0));
}


// ************************************************************************* //
