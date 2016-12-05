/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::hTCModels::noHTC

Description
    Dummy hTC model for 'no hTC'

SourceFiles
    noHTC.C

\*---------------------------------------------------------------------------*/

#ifndef noHTC_H
#define noHTC_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hTCModels
{

/*---------------------------------------------------------------------------*\
                        Class noHTC Declaration
\*---------------------------------------------------------------------------*/

template<class HTempThermoType>
class noHTC
:
    public HTempThermoType
{

    //- Disallow copy construct
    noHTC(const noHTC&);

    //- Disallow default bitwise assignment
    void operator=(const noHTC&);


public:

    //- Runtime type information
    TypeName("noHTC");


    // Constructors

        //- Construct from components
        noHTC
        (
            const word& modelType,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~noHTC();


    // Member Functions

        // Evolution

            //- Correct hTC rate
            virtual void correct();

            //- Fuel consumption rate matrix
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //-  Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

    // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hTCModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "noHTC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
