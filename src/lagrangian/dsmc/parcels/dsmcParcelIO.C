/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "dsmcParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcParcel::dsmcParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    U_(vector::zero),
    RWF_(1.0),
    ERot_(0.0),
    ELevel_(0),
    typeId_(-1),
    newParcel_(0),
    tracked_(nullptr),
    classification_(0),
    stuck_(nullptr),
    vibLevel_(0)
{
    dsmcParcel::TrackedParcel tP = dsmcParcel::TrackedParcel();
    dsmcParcel::StuckParcel sP = dsmcParcel::StuckParcel();
    
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            RWF_ = readScalar(is);
            ERot_ = readScalar(is);
            ELevel_ = readLabel(is);
            typeId_ = readLabel(is);
            newParcel_ = readLabel(is);
            is >> tP;
            classification_ = readLabel(is);
            is >> sP;
            is >> vibLevel_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
                + sizeof(RWF_)
                + sizeof(ERot_)
                + sizeof(ELevel_)
                + sizeof(typeId_)
                + sizeof(newParcel_)
                + sizeof(classification_)
            );
            is >> tP;
            is >> sP;
            is >> vibLevel_;
        }
    }
    
    if (tP.tracked())
    {
        tracked_ = new dsmcParcel::TrackedParcel
            (
                tP.tracked(),
                tP.storePositions(),
                tP.initialTime(),
                tP.initialPosition(),
                tP.meanSquareDisplacementVector(),
                tP.parcelTrajectory()
            );
    }
    
    if (sP.wallTemperature()[0] != 0.0)
    {
        stuck_ = new dsmcParcel::StuckParcel
            (
                sP.wallTemperature(),
                sP.wallVectors()
            );
    }
    
    // Check state of Istream
    is.check
    (
        "Foam::dsmcParcel::dsmcParcel"
        "(const Cloud<dsmcParcel>& cloud, Foam::Istream&), bool"
    );
}


void Foam::dsmcParcel::readFields(Cloud<dsmcParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);
    
    IOField<scalar> RWF
    (
        c.fieldIOobject
        (
            "radialWeight", 
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 1.0)
    );
    c.checkFieldIOobject(c, RWF);

    IOField<scalar> ERot
    (
        c.fieldIOobject
        (
            "ERot", 
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 0.0)
    );
    c.checkFieldIOobject(c, ERot);
    
    IOField<label> ELevel
    (
        c.fieldIOobject
        (
            "ELevel", 
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, ELevel);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newParcel);
    
    IOField<label> classification(c.fieldIOobject("classification", IOobject::MUST_READ));
    c.checkFieldIOobject(c, classification);

    IOField<label> stuckToWall
    (
        c.fieldIOobject
        (
            "stuckToWall", 
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, stuckToWall);
    
    IOField<scalarField> wallTemperature
    (
        c.fieldIOobject
        (
            "wallTemperature", 
            IOobject::READ_IF_PRESENT
        )
    );
    
    if(wallTemperature.size() != c.size())
    {
        wallTemperature.setSize(c.size());
        forAll(wallTemperature, i)
        {
            wallTemperature[i] = scalarField(4, 0.0);
        }
    }
    
    c.checkFieldIOobject(c, wallTemperature);
    
    IOField<vectorField> wallVectors
    (
        c.fieldIOobject
        (
            "wallVectors", 
            IOobject::READ_IF_PRESENT
        )
    );
    
    if(wallVectors.size() != c.size())
    {
        wallVectors.setSize(c.size());
        forAll(wallVectors, i)
        {
            wallVectors[i] = vectorField(4, vector::zero);
        }
    }
    
    c.checkFieldIOobject(c, wallVectors);
    
    IOField<label> isTracked
    (
        c.fieldIOobject
        (
            "isTracked", 
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, isTracked);
    
    IOField<scalar> tracerInitialTime
    (
        c.fieldIOobject
        (
            "tracerInitialTime", 
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 0.0)
    );
    c.checkFieldIOobject(c, tracerInitialTime);
    
    IOField<vector> tracerInitialPosition
    (
        c.fieldIOobject
        (
            "tracerInitialPosition", 
            IOobject::READ_IF_PRESENT
        ),
        vectorField(c.size(), vector::zero)
    );
    c.checkFieldIOobject(c, tracerInitialPosition);
    
    IOField<labelField> vibLevel
    (
        c.fieldIOobject
        (
            "vibLevel", 
            IOobject::READ_IF_PRESENT
        )
    );
    
    if(vibLevel.size() != c.size())
    {
        vibLevel.setSize(c.size());
        forAll(vibLevel, i)
        {
            vibLevel[i].setSize(0);
        }
    }
    
    c.checkFieldIOobject(c, vibLevel);

    label i = 0;
    forAllIter(dsmcCloud, c, iter)
    {
        dsmcParcel& p = iter();

        p.U_ = U[i];
        p.RWF_ = RWF[i];
        p.ERot_ = ERot[i];
        p.ELevel_ = ELevel[i];
        p.typeId_ = typeId[i];
        p.newParcel_ = newParcel[i];
        p.classification_ = classification[i];
        
        if(stuckToWall[i])
        {
            p.setStuck(wallTemperature[i], wallVectors[i]);
        }
        
        if(isTracked[i])
        {
            p.setTracked(true, tracerInitialTime[i], tracerInitialPosition[i]);
        }
        
        p.vibLevel_ = vibLevel[i];
        
        i++;
    }
}


void Foam::dsmcParcel::writeFields(const Cloud<dsmcParcel>& c)
{
    particle::writeFields(c);

    const label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<scalar> ERot(c.fieldIOobject("ERot", IOobject::NO_READ), np);
    IOField<labelField> vibLevel(c.fieldIOobject("vibLevel", IOobject::NO_READ), np);
    IOField<label> ELevel(c.fieldIOobject("ELevel", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::NO_READ), np);
    IOField<label> classification(c.fieldIOobject("classification", IOobject::NO_READ), np);
    
    IOField<label> stuckToWall(c.fieldIOobject("stuckToWall", IOobject::NO_READ), np);
    IOField<scalarField> wallTemperature(c.fieldIOobject("wallTemperature", IOobject::NO_READ), np);
    IOField<vectorField> wallVectors(c.fieldIOobject("wallVectors", IOobject::NO_READ), np);
    
    IOField<label> isTracked(c.fieldIOobject("isTracked", IOobject::NO_READ), np);
    IOField<scalar> tracerInitialTime(c.fieldIOobject("tracerInitialTime", IOobject::NO_READ), np);
    IOField<vector> tracerInitialPosition(c.fieldIOobject("tracerInitialPosition", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(dsmcCloud, c, iter)
    {
        const dsmcParcel& p = iter();

        U[i] = p.U();
        RWF[i] = p.RWF();
        ERot[i] = p.ERot();
        vibLevel[i] = p.vibLevel();
        ELevel[i] = p.ELevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();
        
        stuckToWall[i] = p.isStuck();
        if(stuckToWall[i])
        {
            wallTemperature[i] = p.stuck().wallTemperature();
            wallVectors[i] = p.stuck().wallVectors();
        }
        
        isTracked[i] = p.isTracked();
        if(isTracked[i])
        {
            tracerInitialTime[i] = p.tracked().initialTime();
            tracerInitialPosition[i] = p.tracked().initialPosition();
        }
        
        i++;
    }

    U.write();
    
    if(gMax(RWF) > 1.0)
    {
        //- this is an axisymmetric simulation
        RWF.write();
    }
    
    if(gMax(ERot) > 0.0)
    {
        //- there is at least one molecule
        ERot.write();
    }

    if(gMax(ELevel) > 0)
    {
        //- the electronic mode is activated
        ELevel.write();
    }
    
    typeId.write();
    newParcel.write();
    classification.write();
    
    if(gSum(stuckToWall) > 0)
    {
        //- there is at least one stickingWallPatch with a particle stuck on it
        stuckToWall.write();
        wallTemperature.write(); 
        wallVectors.write();
    }
    
    vibLevel.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const dsmcParcel& p
)
{
    dsmcParcel::TrackedParcel tP = dsmcParcel::TrackedParcel();
    if(p.isTracked())
    {
        tP = p.tracked();
    }
    
    dsmcParcel::StuckParcel sP = dsmcParcel::StuckParcel();
    if(p.isStuck())
    {
        sP = p.stuck();
    }
    
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.RWF()
            << token::SPACE << p.ERot()
            << token::SPACE << p.ELevel()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << tP
            << token::SPACE << p.classification()
            << token::SPACE << sP
            << token::SPACE << p.vibLevel();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            sizeof(p.U())
            + sizeof(p.RWF())
            + sizeof(p.ERot())
            + sizeof(p.ELevel())
            + sizeof(p.typeId())
            + sizeof(p.newParcel())
            + sizeof(p.classification())
        );
        
        os << tP;
        os << sP; 
        os << p.vibLevel();
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::dsmcParcel&)"
    );
    
    return os;
}


// ************************************************************************* //
