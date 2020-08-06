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
    Foam::hTC2Models::noHTC2

Description
    Dummy hTC model for 'no hTC'

SourceFiles
    noHTC2.C

\*---------------------------------------------------------------------------*/

#ifndef noHTC2_H
#define noHTC2_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hTC2Models
{

/*---------------------------------------------------------------------------*\
                        Class noHTC2 Declaration
\*---------------------------------------------------------------------------*/

template<class HTempThermoType>
class noHTC2
:
    public HTempThermoType
{

    //- Disallow copy construct
    noHTC2(const noHTC2&);

    //- Disallow default bitwise assignment
    void operator=(const noHTC2&);


public:

    //- Runtime type information
    TypeName("noHTC2");


    // Constructors

        //- Construct from components
        noHTC2
        (
            const word& modelType,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~noHTC2();


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

            //-  Return chemistry-vibrational coupling source term
            //  in the vibrational equation [kg/m/s3]
            virtual tmp<volScalarField> Scv() const;
            
            //-  Return chemistry-vibrational coupling source term
            //  in the vibrational equation for specie i [kg/m/s3]
            virtual tmp<volScalarField> Scv(const label i) const;

            //- Return energy removed from free-electrons during impact
            //  ionisation reactions [kg/m/s3]
            virtual tmp<volScalarField> Seiir() const;
            
            //- Return energy removed from free-electrons during impact
            //  ionisation reactions for specie i [kg/m/s3]
            virtual tmp<volScalarField> Seiir(const label i) const;

    // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hTC2Models
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "noHTC2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
