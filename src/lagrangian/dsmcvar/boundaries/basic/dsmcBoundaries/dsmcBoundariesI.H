/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline List< autoPtr<dsmcPatchBoundary> >& dsmcBoundaries::patchBoundaryModels()
{
    return patchBoundaryModels_;
}

inline const List< autoPtr<dsmcPatchBoundary> >& dsmcBoundaries::patchBoundaryModels() const
{
    return patchBoundaryModels_;
}

inline const List<word>& dsmcBoundaries::patchBoundaryNames() const
{
    return patchBoundaryNames_;
}

inline const List<label>& dsmcBoundaries::patchBoundaryIds() const
{
    return patchBoundaryIds_;
}

inline const List<label>& dsmcBoundaries::patchToModelIds() const
{
    return patchToModelId_;
}




inline List< autoPtr<dsmcCyclicBoundary> >& dsmcBoundaries::cyclicBoundaryModels()
{
    return cyclicBoundaryModels_;
}

inline const List< autoPtr<dsmcCyclicBoundary> >& dsmcBoundaries::cyclicBoundaryModels() const
{
    return cyclicBoundaryModels_;
}

inline const List<word>& dsmcBoundaries::cyclicBoundaryNames() const
{
    return cyclicBoundaryNames_;
}

inline const List<label>& dsmcBoundaries::cyclicBoundaryIds() const
{
    return cyclicBoundaryIds_;
}

inline const List<label>& dsmcBoundaries::cyclicBoundaryToModelIds() const
{
    return cyclicBoundaryToModelId_;
}






inline List< autoPtr<dsmcGeneralBoundary> >& dsmcBoundaries::generalBoundaryModels()
{
    return generalBoundaryModels_;
}

inline const List< autoPtr<dsmcGeneralBoundary> >& dsmcBoundaries::generalBoundaryModels() const
{
    return generalBoundaryModels_;
}

inline const List<word>& dsmcBoundaries::generalBoundaryNames() const
{
    return generalBoundaryNames_;
}

inline const List<label>& dsmcBoundaries::generalBoundaryIds() const
{
    return generalBoundaryIds_;
}








// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
