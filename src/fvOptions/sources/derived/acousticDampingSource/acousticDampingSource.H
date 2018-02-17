/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fv::acousticDampingSource

Group
    grpFvOptionsSources

Description
    Acoustic damping source

    \heading Source usage

    Example usage:
    \verbatim
    acousticDampingSourceCoeffs
    {
    }
    \endverbatim

SourceFiles
    acousticDampingSource.C

\*---------------------------------------------------------------------------*/

#ifndef acousticDampingSource_H
#define acousticDampingSource_H

#include "cellSetOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                    Class acousticDampingSource Declaration
\*---------------------------------------------------------------------------*/

class acousticDampingSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Frequency [Hz]
        dimensionedScalar frequency_;

        //- Blending factor []
        volScalarField blendFactor_;

        //- Name of reference velocity field
        word URefName_;

        // Sphere centre location or damping
        point x0_;

        // Inner radius at which to start damping
        scalar r1_;

        // Outer radius beyond which damping is applied
        scalar r2_;

        //- Stencil width, default = 20
        label w_;


    // Protected Member Functions

        //- Helper function to set the blending factor
        void setBlendingFactor();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        acousticDampingSource(const acousticDampingSource&);

        //- Disallow default bitwise assignment
        void operator=(const acousticDampingSource&) = delete;


public:

    //- Runtime type information
    TypeName("acousticDampingSource");


    // Constructors

        //- Construct from components
        acousticDampingSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~acousticDampingSource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
