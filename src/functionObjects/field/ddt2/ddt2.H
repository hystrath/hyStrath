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
    Foam::functionObjects::ddt2

Group
    grpFieldFunctionObjects

Description
    This function object calculates the magnitude squared
    of d(scalarField)/dt.

    The result can be used further for determining variance or RMS values
    (for example).

Usage
    Example of function object specification:
    \verbatim
    dpdt2
    {
        type        ddt2;
        libs        ("libfieldFunctionObjects.so");
        fields      (p);
        result      d@@dt2;
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required  | Default value
        type     | type name: ddt2            | yes       |
        fields   | Name of fields to process  | yes       |
        result   | Name of results            | no        | magSqr(ddt(@@))
        log      | Log to standard output     | no        | yes
        mag      | Use 'mag' instead of 'magSqr' | no     | false
    \endtable

    Note that the optional 'mag' entry cannot be changed during the simulation
    since it alters the dimensions of the output field.

    A list of fields can contain exact names or regular expressions.
    The token '\@\@' in the result name is replaced by the name of the source
    field.

    The function object will skip over fields that appear to have
    already been processed (ie, their names are similar to the output names).

SourceFiles
    ddt2.C
    ddt2Templates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ddt2_H
#define functionObjects_ddt2_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "OFstream.H"
#include "wordReList.H"
#include "regExp.H"
#include "HashSet.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class ddt2 Declaration
\*---------------------------------------------------------------------------*/

class ddt2
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of fields to process
        wordReList selectFields_;

        //- Formatting for the result fields.
        word resultName_;

        //- Avoid processing the same field twice.
        mutable regExp blacklist_;

        //- Hashed names of result fields.
        wordHashSet results_;

        //- Use 'mag' instead of 'magSqr'.
        //  Cannot be adjusted during the simulation since it alters the
        //  dimensions of the output field.
        const Switch mag_;


    // Private Member Functions

        //- Check that the word contains the appropriate substitution token(s).
        static bool checkFormatName(const word&);


        //- Accept unless field name appears to have already been processed
        bool accept(const word& fieldName) const;

        //- Apply for the volume field type
        template<class FieldType>
        int apply(const word& inputName, int& state);

        //- Process by trying to apply for various volume field types.
        int process(const word& inputName);


        //- Disallow default bitwise copy construct
        ddt2(const ddt2&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ddt2&) = delete;


public:

    //- Runtime type information
    TypeName("ddt2");


    // Constructors

        //- Construct from Time and dictionary
        ddt2
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ddt2();


    // Member Functions

        //- Read the ddt2 specification
        virtual bool read(const dictionary&);

        //- Calculate the ddt2 fields
        virtual bool execute();

        //- Write the ddt fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ddt2Templates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
