/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "fvMesh.H"

template<class Type>
Type Foam::functionObjects::extractEulerianParticles::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
    const label localFacei,
    const label globalFacei
) const
{
    if (mesh_.isInternalFace(globalFacei))
    {
        return field[globalFacei];
    }
    else
    {
        label patchi = patchIDs_[localFacei];
        label pFacei = patchFaceIDs_[localFacei];
        if (patchi != 0)
        {
            return field.boundaryField()[patchi][pFacei];
        }
        else
        {
            return pTraits<Type>::zero;
        }
    }
}


// ************************************************************************* //
