/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
