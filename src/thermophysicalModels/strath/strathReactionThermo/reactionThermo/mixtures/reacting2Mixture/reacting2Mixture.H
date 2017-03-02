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
    Foam::reacting2Mixture

Description
    Foam::reacting2Mixture

SourceFiles
    reacting2Mixture.C

\*---------------------------------------------------------------------------*/

#ifndef reacting2Mixture_H
#define reacting2Mixture_H

#include "speciesTable.H"
#include "chemistry2Reader.H"
#include "multi2ComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class reacting2Mixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class reacting2Mixture
:
    public speciesTable,
    public autoPtr<chemistry2Reader<ThermoType> >,
    public multi2ComponentMixture<ThermoType>,
    public PtrList<Reaction2<ThermoType> >
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        reacting2Mixture(const reacting2Mixture&);

        //- Disallow default bitwise assignment
        void operator=(const reacting2Mixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        reacting2Mixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~reacting2Mixture()
    {}


    // Member functions

        //- Read dictionary
        void read(const dictionary&);

        label size() const
        {
            return PtrList<Reaction2<ThermoType> >::size();
        }

        Reaction2<ThermoType>& operator [] (const label i)
        {
            return PtrList<Reaction2<ThermoType> >::operator[](i);
        }

        const Reaction2<ThermoType>& operator [] (const label i) const
        {
            return PtrList<Reaction2<ThermoType> >::operator[](i);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reacting2Mixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
