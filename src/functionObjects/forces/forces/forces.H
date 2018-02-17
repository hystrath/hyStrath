/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::functionObjects::forces

Group
    grpForcesFunctionObjects

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches, and the resistance
    from porous zones.

    Forces and moments are calculated, with optional co-ordinate system and
    writing of binned data, where force and moment contributions are collected
    into a user-defined number of bins that span the input geometries for a
    user-defined direction vector.

    Data is written into multiple files in the
    postProcessing/\<functionObjectName\> directory:
    - force.dat          : forces in global Cartesian co-ordinate system
    - moment.dat         : moments in global Cartesian co-ordinate system
    - forceBin.dat       : force bins in global Cartesian co-ordinate system
    - momentBin.dat      : moment bins in global Cartesian co-ordinate system
    - localForce.dat     : forces in local co-ordinate system
    - localMoment.dat    : moments in local co-ordinate system
    - localForceBin.dat  : force bins in local co-ordinate system
    - localMomentBin.dat : moment bins in local co-ordinate system

Usage
    Example of function object specification:
    \verbatim
    forces1
    {
        type        forces;
        libs        ("libstrathforces.so");
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: forces       | yes         |
        log          | Write force data to standard output | no | no
        writeFields  | Write the force and moment fields | no | no
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        CofR         | Centre of rotation (see below) | no   |
        porosity     | flag to include porosity contributions | no | no
        directForceDensity | Force density supplied directly (see below)|no|no
        fD           | Name of force density field (see below) | no | fD
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with incresing distance | yes |
    \endtable

Note
  - For incompressible cases, set \c rho to \c rhoInf.  You will then be
    required to provide a \c rhoInf value corresponding to the free-stream
    constant density.
  - If the force density is supplied directly, set the \c directForceDensity
    flag to 'yes', and supply the force density field using the \c
    fDName entry
  - The centre of rotation (CofR) for moment calculations can either be
    specified by an \c CofR entry, or be taken from origin of the local
    coordinate system.  For example,
    \verbatim
        CofR        (0 0 0);
    \endverbatim
    or
    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            e3      (0 0 1);
            e1      (1 0 0);
        }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forces_H
#define functionObjects_forces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class forces
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;

        // File streams

            //- Forces
            autoPtr<OFstream> forceFilePtr_;

            //- Moments
            autoPtr<OFstream> momentFilePtr_;

            //- Force bins
            autoPtr<OFstream> forceBinFilePtr_;

            //- Moment bins
            autoPtr<OFstream> momentBinFilePtr_;

            //- Local force
            autoPtr<OFstream> localForceFilePtr_;

            //- Local moment
            autoPtr<OFstream> localMomentFilePtr_;

            //- Local force bins
            autoPtr<OFstream> localForceBinFilePtr_;

            //- Local moment bins
            autoPtr<OFstream> localMomentBinFilePtr_;


        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- Is the force density being supplied directly?
            Switch directForceDensity_;

            //- The name of the force density (fD) field
            word fDName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Coordinate system used when evaluting forces/moments
            coordinateSystem coordSys_;

            //- Flag to indicate whether we are using a local co-ordinate sys
            bool localSystem_;

            //- Flag to include porosity effects
            bool porosity_;


            // Bin information

                //- Number of bins
                label nBin_;

                //- Direction used to determine bin orientation
                vector binDir_;

                //- Distance between bin divisions
                scalar binDx_;

                //- Minimum bin bounds
                scalar binMin_;

                //- Bin positions along binDir
                List<point> binPoints_;

                //- Should bin data be cumulative?
                bool binCumulative_;


            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;


    // Protected Member Functions

        //- Create a field name
        word fieldName(const word& name) const;

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Write header for binned data
        void writeBinHeader(const word& header, Ostream& os) const;

        //- Initialise the fields
        void initialise();

        //- Initialise the collection bins
        void initialiseBins();

        //- Reset the fields prior to accumulation of force/moments
        void resetFields();

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Accumulate bin data
        void applyBins
        (
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP,
            const vectorField& d
        );

        //- Add patch contributions to force and moment fields
        void addToFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Add cell contributions to force and moment fields
        void addToFields
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Helper function to write integrated forces and moments
        void writeIntegratedForceMoment
        (
            const string& descriptor,
            const vectorField& fm0,
            const vectorField& fm1,
            const vectorField& fm2,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write force data
        void writeForces();

        //- Helper function to write binned forces and moments
        void writeBinnedForceMoment
        (
            const List<Field<vector>>& fm,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write binned data
        void writeBins();

        //- Disallow default bitwise copy construct
        forces(const forces&);

        //- Disallow default bitwise assignment
        void operator=(const forces&);


public:

    //- Runtime type information
    TypeName("forces");


    // Constructors

        //- Construct from Time and dictionary
        forces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        forces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~forces();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Calculate the forces and moments
        virtual void calcForcesMoment();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
