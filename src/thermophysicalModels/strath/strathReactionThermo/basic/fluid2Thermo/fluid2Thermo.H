/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::fluid2Thermo

Description
    Fundamental fluid thermodynamic properties

SourceFiles
    fluid2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluid2Thermo_H
#define fluid2Thermo_H

#include "basic2Thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluid2Thermo Declaration
\*---------------------------------------------------------------------------*/

class fluid2Thermo
:
    public basic2Thermo
{

protected:

    // Protected data

        //- Transport dictionary
        IOdictionary transportDictionary_;
        

public:

    //- Runtime type information
    TypeName("fluid2Thermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluid2Thermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Constructors

        //- Construct from mesh and phase name
        fluid2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        fluid2Thermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<fluid2Thermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~fluid2Thermo();


    // Member functions

        static const word transportDictName;
        
        IOdictionary& transportDictionary()
        {
            return transportDictionary_;
        }
        
        // Access to thermodynamic state variables

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;
            
            //- Return non-const access to the local compressibility 
            //  field [s^2/m^2]
            virtual volScalarField& getPsi() = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual const scalarField& mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
