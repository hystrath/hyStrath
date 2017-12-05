/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::valueAverage::calc
(
    const word& fieldName,
    const word& meanName,
    const scalar alpha,
    const scalar beta,
    bool& processed
)
{
    const word valueType = objectResultType(functionObjectName_, fieldName);

    if (pTraits<Type>::typeName != valueType)
    {
        return;
    }

    Type currentValue = getObjectResult<Type>(functionObjectName_, fieldName);

    Type meanValue = getResult<Type>(meanName);
    meanValue = alpha*meanValue + beta*currentValue;

    setResult(meanName, meanValue);

    file() << tab << meanValue;

    Log<< "    " << meanName << ": " << meanValue << nl;

    processed = true;
}


// ************************************************************************* //
