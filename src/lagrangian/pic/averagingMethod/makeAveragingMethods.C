/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "fvCFD.H"

#include "polyMeshTetDecomposition.H"

#include "Basic.H"
#include "Linear.H"


// Scalar interpolation
defineNamedTemplateTypeNameAndDebug(Foam::AveragingMethod<Foam::scalar>, 0);
namespace Foam
{
    defineTemplateRunTimeSelectionTable
    (
        AveragingMethod<Foam::scalar>,
        dictionary
    );
}

// Vector interpolation
defineNamedTemplateTypeNameAndDebug(Foam::AveragingMethod<Foam::vector>, 0);
namespace Foam
{
    defineTemplateRunTimeSelectionTable
    (
        Foam::AveragingMethod<Foam::vector>,
        dictionary
    );
}


// Basic interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Basic<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Basic<Foam::scalar>>
    addBasicscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Basic<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Basic<Foam::vector>>
    addBasicvectorConstructorToTable_;
    
    
// Linear interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Linear<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Linear<Foam::scalar>>
    addLinearscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Linear<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Linear<Foam::vector>>
    addLinearvectorConstructorToTable_;    
    

/*namespace Foam
{
    namespace AveragingMethods
    {
        // Basic interpolation
        defineNamedTemplateTypeNameAndDebug(Basic<scalar>, 0);
        AveragingMethod<scalar>::
            adddictionaryConstructorToTable<Basic<scalar> >
            addBasicscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Basic<vector>, 0);
        AveragingMethod<vector>::
            adddictionaryConstructorToTable<Basic<vector> >
            addBasicvectorConstructorToTable_;

        // Dual interpolation // Chris note: not working 
        //defineNamedTemplateTypeNameAndDebug(Dual<scalar>, 0);
        //AveragingMethod<scalar>::
        //    adddictionaryConstructorToTable<Dual<scalar> >
        //    addDualscalarConstructorToTable_;

        //defineNamedTemplateTypeNameAndDebug(Dual<vector>, 0);
        //AveragingMethod<vector>::
        //    adddictionaryConstructorToTable<Dual<vector> >
        //    addDualvectorConstructorToTable_;
        
        // Linear interpolation
        defineNamedTemplateTypeNameAndDebug(Linear<scalar>, 0);
        AveragingMethod<scalar>::
            adddictionaryConstructorToTable<Linear<scalar> >
            addLinearscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Linear<vector>, 0);
        AveragingMethod<vector>::
            adddictionaryConstructorToTable<Linear<vector> >
            addLinearvectorConstructorToTable_;
    
}*/


// ************************************************************************* //
