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

Class


Description

\*---------------------------------------------------------------------------*/

template<class ParticleType>
inline const Foam::labelListList& Foam::cellInteractions<ParticleType>::inverseDIL() const
{
    return inverseDIL_;
}

template<class ParticleType>
inline const Foam::labelListList& Foam::cellInteractions<ParticleType>::dil() const
{
    return dil_;
}

template<class ParticleType>
inline const Foam::labelListList& Foam::cellInteractions<ParticleType>::fil() const
{
    return fil_;
}

template<class ParticleType>
inline const Foam::List<Foam::referredCell>& Foam::cellInteractions<ParticleType>::refCells() const
{
    return refCells_;
}

template<class ParticleType>
inline Foam::List<Foam::referredCell>& Foam::cellInteractions<ParticleType>::refCells()
{
    return refCells_;
}

template<class ParticleType>
inline const Foam::List<Foam::DynamicList<ParticleType*> >&
    Foam::cellInteractions<ParticleType>::refCellsParticles() const
{
    return refCellsParticles_;
}

template<class ParticleType>
inline Foam::List<Foam::DynamicList<ParticleType*> >&
    Foam::cellInteractions<ParticleType>::refCellsParticles()
{
    return refCellsParticles_;
}

template<class ParticleType>
inline const Foam::List<Foam::DynamicList<ParticleType*> >&
    Foam::cellInteractions<ParticleType>::ripl() const
{
    return ripl_;
}

template<class ParticleType>
inline const Foam::labelListList& Foam::cellInteractions<ParticleType>::sourceCellToRefs() const
{
    return sourceCellToRefs_;
}

template<class ParticleType>
inline const Foam::labelListList& Foam::cellInteractions<ParticleType>::inverseFRIL() const
{
    return inverseFRIL_;
}

template<class ParticleType>
inline Foam::IDLList<ParticleType>& Foam::cellInteractions<ParticleType>::referredCloud()
{
    return referredCloud_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
