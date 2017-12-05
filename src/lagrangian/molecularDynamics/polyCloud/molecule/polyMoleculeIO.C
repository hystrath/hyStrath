/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "polyMolecule.H"
#include "IOstreams.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMolecule::polyMolecule
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    Q_(tensor::zero),
    v_(vector::zero),
    a_(vector::zero),
    pi_(vector::zero),
    tau_(vector::zero),
    specialPosition_(vector::zero),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(0),
    id_(0),
    R_(GREAT),
    frac_(1.0),
    trackingNumber_(-1),
    siteForces_(0),
    sitePositions_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> Q_;
            is  >> v_;
            is  >> a_;
            is  >> pi_;
            is  >> tau_;
            is  >> specialPosition_;
            potentialEnergy_ = readScalar(is);
            is  >> rf_;
            special_ = readLabel(is);
            id_ = readLabel(is);
            is >> R_;
            is >> frac_;
            is >> trackingNumber_;
            is  >> siteForces_;
            is  >> sitePositions_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&Q_),
                sizeof(Q_)
              + sizeof(v_)
              + sizeof(a_)
              + sizeof(pi_)
              + sizeof(tau_)
              + sizeof(specialPosition_)
              + sizeof(potentialEnergy_)
              + sizeof(rf_)
              + sizeof(special_)
              + sizeof(id_)
              + sizeof(R_)
              + sizeof(frac_)
              + sizeof(trackingNumber_)
            );

            is >> siteForces_ >> sitePositions_;
        }
    }

    // Check state of Istream
    is.check
    (
        "Foam::polyMolecule::polyMolecule"
        "(const Cloud<polyMolecule>& cloud, Foam::Istream&), bool"
    );
}


void Foam::polyMolecule::readFields(Cloud<polyMolecule>& mC)
{
    if (!mC.size())
    {
        return;
    }

    particle::readFields(mC);
    
    Info << "Reading fields" << endl;

    IOField<tensor> Q(mC.fieldIOobject("Q", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, Q);

    // MB: I have removed reading/writing these fields because no one ever uses them
    //    and they by doing this we are saving a lot on hard disk space
    
//     IOField<tensor> rf(mC.fieldIOobject("rf", IOobject::MUST_READ));
//     mC.checkFieldIOobject(mC, rf);

    IOField<vector> v(mC.fieldIOobject("v", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, v);

//     IOField<vector> a(mC.fieldIOobject("a", IOobject::MUST_READ));
//     mC.checkFieldIOobject(mC, a);

    IOField<vector> pi(mC.fieldIOobject("pi", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, pi);

//     IOField<vector> tau(mC.fieldIOobject("tau", IOobject::MUST_READ));
//     mC.checkFieldIOobject(mC, tau);

    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::MUST_READ)
    );
    
    mC.checkFieldIOobject(mC, specialPosition);

    IOField<label> special(mC.fieldIOobject("special", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, special);

    IOField<label> id(mC.fieldIOobject("id", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, id);
    
    IOField<label> trackingNumber(mC.fieldIOobject("trackingNumber", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, trackingNumber);    

    label i = 0;
    forAllIter(polyMoleculeCloud, mC, iter)
    {
        polyMolecule& mol = iter();

        mol.Q_ = Q[i];
//         mol.rf_ = rf[i];
        mol.v_ = v[i];
//         mol.a_ = a[i];
        mol.pi_ = pi[i];
//         mol.tau_ = tau[i];
        mol.specialPosition_ = specialPosition[i];
        mol.special_ = special[i];
        mol.id_ = id[i];
        mol.trackingNumber_ = trackingNumber[i];        
        i++;
    }
}


void Foam::polyMolecule::writeFields(const Cloud<polyMolecule>& mC)
{
    particle::writeFields(mC);
    
    label np = mC.size();

    IOField<tensor> Q(mC.fieldIOobject("Q", IOobject::NO_READ), np);
//     IOField<tensor> rf(mC.fieldIOobject("rf", IOobject::NO_READ), np);
    IOField<vector> v(mC.fieldIOobject("v", IOobject::NO_READ), np);
//     IOField<vector> a(mC.fieldIOobject("a", IOobject::NO_READ), np);
    IOField<vector> pi(mC.fieldIOobject("pi", IOobject::NO_READ), np);
//     IOField<vector> tau(mC.fieldIOobject("tau", IOobject::NO_READ), np);
    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::NO_READ),
        np
    );
    IOField<label> special(mC.fieldIOobject("special", IOobject::NO_READ), np);
    IOField<label> id(mC.fieldIOobject("id", IOobject::NO_READ), np);
    IOField<label> trackingNumber(mC.fieldIOobject("trackingNumber", IOobject::NO_READ), np);
    
    // Post processing fields

    // MB: I have removed reading/writing these fields because no one ever uses them
    //    and they by doing this we are saving a lot on hard disk space
    
//     IOField<vector> piGlobal
//     (
//         mC.fieldIOobject("piGlobal", IOobject::NO_READ),
//         np
//     );
// 
//     IOField<vector> tauGlobal
//     (
//         mC.fieldIOobject("tauGlobal", IOobject::NO_READ),
//         np
//     );
// 
//     IOField<vector> orientation1
//     (
//         mC.fieldIOobject("orientation1", IOobject::NO_READ),
//         np
//     );
// 
//     IOField<vector> orientation2
//     (
//         mC.fieldIOobject("orientation2", IOobject::NO_READ),
//         np
//     );
// 
//     IOField<vector> orientation3
//     (
//         mC.fieldIOobject("orientation3", IOobject::NO_READ),
//         np
//     );

    label i = 0;
    forAllConstIter(polyMoleculeCloud, mC, iter)
    {
        const polyMolecule& mol = iter();

//         rf[i] = mol.rf_;
        Q[i] = mol.Q_;
        v[i] = mol.v_;
//         a[i] = mol.a_;
        pi[i] = mol.pi_;
//         tau[i] = mol.tau_;
        specialPosition[i] = mol.specialPosition_;
        special[i] = mol.special_;
        id[i] = mol.id_;
        trackingNumber[i] = mol.trackingNumber_;
        
//         piGlobal[i] = mol.Q_ & mol.pi_;
//         tauGlobal[i] = mol.Q_ & mol.tau_;

//         orientation1[i] = mol.Q_ & vector(1,0,0);
//         orientation2[i] = mol.Q_ & vector(0,1,0);
//         orientation3[i] = mol.Q_ & vector(0,0,1);

        i++;
    }

//     rf.write();
    Q.write();
    v.write();
//     a.write();
    pi.write();
//     tau.write();
    specialPosition.write();
    special.write();
    id.write();
    trackingNumber.write();
    
//     piGlobal.write();
//     tauGlobal.write();

//     orientation1.write();
//     orientation2.write();
//     orientation3.write();

    Info<< "writeFields " << mC.name() << endl;
    
    if (isA<polyMoleculeCloud>(mC))
    {
        const polyMoleculeCloud& m = dynamic_cast<const polyMoleculeCloud&>(mC);

        m.writeXYZ
        (
            m.mesh().time().timePath()/cloud::prefix/"polyMoleculeCloud.xmol"
        );
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyMolecule& mol)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << token::SPACE << static_cast<const particle&>(mol)
            << token::SPACE << mol.face()
            << token::SPACE << mol.stepFraction()
            << token::SPACE << mol.Q_
            << token::SPACE << mol.v_
            << token::SPACE << mol.a_
            << token::SPACE << mol.pi_
            << token::SPACE << mol.tau_
            << token::SPACE << mol.specialPosition_
            << token::SPACE << mol.potentialEnergy_
            << token::SPACE << mol.rf_
            << token::SPACE << mol.special_
            << token::SPACE << mol.id_
            << token::SPACE << mol.R_
            << token::SPACE << mol.frac_
            << token::SPACE << mol.trackingNumber_
            << token::SPACE << mol.siteForces_
            << token::SPACE << mol.sitePositions_;
    }
    else
    {
        os  << static_cast<const particle&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.Q_),
            sizeof(mol.Q_)
          + sizeof(mol.v_)
          + sizeof(mol.a_)
          + sizeof(mol.pi_)
          + sizeof(mol.tau_)
          + sizeof(mol.specialPosition_)
          + sizeof(mol.potentialEnergy_)
          + sizeof(mol.rf_)
          + sizeof(mol.special_)
		  + sizeof(mol.id_)
          + sizeof(mol.R_)
          + sizeof(mol.frac_)
          + sizeof(mol.trackingNumber_)
        );
        os << mol.siteForces_ << mol.sitePositions_;
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::polyMolecule&)"
    );

    return os;
}


// ************************************************************************* //
