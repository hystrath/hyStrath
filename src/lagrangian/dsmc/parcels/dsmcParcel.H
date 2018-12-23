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

Class
    Foam::dsmcParcel

Description
    DSMC parcel class

SourceFiles
    dsmcParcelI.H
    dsmcParcel.C
    dsmcParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcParcel_H
#define dsmcParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dsmcCloud;

// Forward declaration of friend functions
/*Ostream& operator<<
(
    Ostream&,
    const dsmcParcel&
)*/


/*---------------------------------------------------------------------------*\
                         Class dsmcParcel Declaration
\*---------------------------------------------------------------------------*/

class dsmcParcel
:
      public particle
{
public:

    //- Class to hold DSMC particle constant properties
    class constantProperties
    {
        // Private data

            //- Particle type deduced from its other properties
            //  0: electron, 10: atom, 11: ionised atom, 
            //  20: neutral diatomic molecule, 21: ionised diatomic molecule
            //  30: neutral polyatomic molecule, 31: ionised polyatomic molecule
            label type_;
            
            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Rotational degrees of freedom
            scalar rotationalDegreesOfFreedom_;
        
            //- Number of vibrational modes
            label nVibrationalModes_;
            
            //- Degeneracy of each vibrational mode
            labelList vibrationalDegeneracyList_;

            //- Viscosity index
            scalar omega_;
            
            //- VSS scattering parameter
            scalar alpha_;
        
            //- Characteristic vibrational temperature
            scalarList thetaV_;

            //- Characteristic dissociation temperature
            scalarList thetaD_;
        
            //- Reference vibrational relaxation number
            scalarList Zref_;
            
            //- Reference temperature for reference vibrational relaxation 
            //  number
            scalarList TrefZv_;
            
            //- Characteristic dissociation quantum level
            labelList charDissQuantumLevel_;
            
            //- Ionisation temperature
            scalar ionisationTemperature_;
            
            //- Charge state (0 for neutrals, 1 for ions, -1 for electrons)
            label charge_;
            
            //- Number of electronic energy/degeneracy levels
            label nElectronicLevels_;
            
            //- Scalar list of electronic energy levels [J]
            scalarList electronicEnergyList_;     
        
            //- Degeneracy of each electronic level
            labelList electronicDegeneracyList_;

    public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member functions

            //- Return const access to the type of particle
            inline label type() const;
            
            //- Return const access to the particle mass
            inline scalar mass() const;

            //- Return const access to the particle diameter
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the rotational degrees of freedom
            inline scalar rotationalDegreesOfFreedom() const;
            
            //- Return the number of vibrational modes
            inline label nVibrationalModes() const;
            
            //- Return the degeneracy for the vibrational mode m
            inline label vibrationalDegeneracy_m(const label mode) const;
            
            //- Return the degeneracy of each vibrational mode
            inline labelList vibrationalDegeneracyList() const;
            
            //- Return the viscosity index
            inline scalar omega() const;
            
            //- Return the VSS scattering parameter
            inline scalar alpha() const;
    
            //- Return the characteristic vibrational temperature for a given
            //  vibrational energy mode
            inline scalar thetaV_m(const label mode) const;
            
            //- Return the characteristic vibrational temperature
            inline scalarList thetaV() const;    
            
            //- Return the characteristic dissociation temperature for a given
            //  vibrational energy mode
            inline scalar thetaD_m(const label mode) const;
            
            //- Return the characteristic dissociation temperature
            inline scalarList thetaD() const;
    
            //- Return the characteristic vibrational relaxation number for a 
            //  given vibrational energy mode
            inline scalar Zref_m(const label mode) const;
            
            //- Return the characteristic vibrational relaxation number
            inline scalarList Zref() const;
                      
            //- Return the reference temperature for the characteristic 
            //  vibrational relaxation number for a given vibrational energy 
            //  mode
            inline scalar TrefZv_m(const label mode) const;
            
            //- Return the reference temperature for the characteristic 
            //  vibrational relaxation number
            inline scalarList TrefZv() const;
            
            //- Return the characteristic dissociation level for a given
            //  vibrational energy mode
            inline label charDissQuantumLevel_m(const label mode) const;
            
            //- Return the characteristic dissociation level
            inline labelList charDissQuantumLevel() const;
            
            //- Return the vibrational energy of a given vibrational energy mode
            inline scalar eVib_m(const label mode, const label vibLevel) const;
            
            //- Return the total vibrational energy given the vibrational level 
            //  for each mode
            inline scalar eVib_tot(const labelList& vibLevel) const;
            
            //- Return the vibrational energy for each vibrational energy mode
            //  given the vibrational level for each mode
            inline scalarList eVib(const labelList& vibLevel) const;
            
            //- Return the ionisation temperature
            inline scalar ionisationTemperature() const;
            
            //- Return the particle's charge state
            inline label charge() const;
            
            //- Return the number of electronic energy levels
            inline label nElectronicLevels() const;
            
            //- Return the energy (J) associated to each electronic level
            inline scalarList electronicEnergyList() const;
            
            //- Return the degeneracy of each electronic level
            inline labelList electronicDegeneracyList() const;
    };


    //- Class used to store info about parcels stuck on a sticking wall
    class StuckParcel
    {
        protected:
        
            // Protected data
            
              //- Store temperature of wall particle is stuck to
              //  [0]: local temperature, [1]: patch id, [2]: local face,
              //  [3]: pre-interaction energy
              scalarField wallTemperature_;
              
              //- Store vectors of wall particle is stuck to
              //  [0 and 1]: first and second wall unit tangential vectors, 
              //  [2]: wall unit normal vector, [3]: pre-interaction momentum
              vectorField wallVectors_;
            
            
            // Protected member functions
            
              //- Construct and return a clone
              StuckParcel* clone() const
              {
                  return new StuckParcel(*this);
              }
            
            
        public:
        
            // Constructors

                //- Construct from components
                StuckParcel()
                :
                    wallTemperature_(scalarField(4, 0.0)),
                    wallVectors_(vectorField(4, vector::zero))
                {}
                
                StuckParcel
                (
                    const scalarField& wallTemperature,
                    const vectorField& wallVector
                )
                :
                    wallTemperature_(wallTemperature),
                    wallVectors_(wallVector)
                {}
                
                //- Copy constructor
                StuckParcel
                (
                    const StuckParcel& sP
                )
                :
                    wallTemperature_(sP.wallTemperature_),
                    wallVectors_(sP.wallVectors_)
                {}
                
                
            // Destructor
                
                ~StuckParcel(){}
                
                
            // Member functions
            
            
            
            // Access
            
                //- Return const access to wallTemperature
                inline const scalarField& wallTemperature() const
                {
                    return wallTemperature_;
                }
                
                //- Return const access to wallVectors
                inline const vectorField& wallVectors() const
                {
                    return wallVectors_;
                }
            
            
            // Edit
            
                //- Return access to wallTemperature
                inline scalarField& wallTemperature()
                {
                    return wallTemperature_;
                }
                
                //- Return access to wallVectors
                inline vectorField& wallVectors()
                {
                    return wallVectors_;
                }
                
                
            // Operators
                
                //- Assignment operator
                StuckParcel& operator=(const StuckParcel& sP)
                {
                     wallTemperature_ = sP.wallTemperature_;
                     wallVectors_ = sP.wallVectors_;
                     
                     return *this;
                }
                
                //- Istream operator
                friend Istream& operator>>         
                (
                    Istream& is,
                    StuckParcel& sP
                )
                {
                    if (is.format() == IOstream::ASCII)
                    {
                        is >> sP.wallTemperature_;
                        is >> sP.wallVectors_;
                    }
                    else
                    {
                        is >> sP.wallTemperature_;
                        is >> sP.wallVectors_;
                    }
                    
                    // Check state of Istream
                    is.check
                    (
                        "Foam::Istream& Foam::operator>>"
                        "(Foam::Istream&, Foam::dsmcParcel::TrackedParcel&)"
                    );
                    
                    return is;
                }
                
                //- Ostream operator
                friend Ostream& operator<<         
                (
                    Ostream& os,
                    const StuckParcel& sP
                )
                {
                    if (os.format() == IOstream::ASCII)
                    {
                        os << sP.wallTemperature_ << token::SPACE
                           << sP.wallVectors_;
                    }
                    else
                    {
                        os << sP.wallTemperature_
                           << sP.wallVectors_;
                    }
                       
                    // Check state of Ostream
                    os.check
                    (
                        "Foam::Ostream& Foam::operator<<"
                        "(Foam::Ostream&, const Foam::dsmcParcel::StuckParcel&)"
                    );
                    
                    return os;
                }

    };        
            
    //- Class used to track randomly-selected parcels
    class TrackedParcel
    {
        public:
        
            // Static data
        
              //- Number of tracked parcels having exited the domain
              static label nDELETED;
              

        protected:
        
            // Protected data
            
              //- Is the parcel being tracked
              bool tracked_;
              
              //- Patch iD of the boundary where the parcel is being added
              label inPatchId_;
              
              //- Record the parcel's trajectory
              bool storePositions_;
              
              //- Time at which the parcel was selected as a tracer
              scalar t0_;
              
              //- Parcel's initial position
              vector initialPosition_;
              
              //- Parcel's current position
              vector currentPosition_;
              
              //- Distance travelled
              vector distanceTravelled_;
              
              //- Mean squared displacement
              vector msd_;
              
              //- A data structure holding the position of the tracked parcel
              //  versus time
              //DynamicList<FixedList<scalar, 4>> parcelTrajectory_;
            
            
            // Protected member functions
            
              //- Construct and return a clone
              TrackedParcel* clone() const
              {
                  return new TrackedParcel(*this);
              }
            
            
            
        public:
        
            // Constructors

                //- Construct from components
                TrackedParcel()
                :
                    tracked_(false),
                    inPatchId_(-1),
                    storePositions_(false),
                    t0_(0.0),
                    initialPosition_(vector::zero),
                    currentPosition_(vector::zero),
                    distanceTravelled_(vector::zero),
                    msd_(vector::zero)//,
                    //parcelTrajectory_()
                {
                    //Info << "Construct from components" << endl;
                }
                
                TrackedParcel
                (
                    const label& inPatchId,
                    const scalar& initialTime,
                    const vector& initialPosition
                )
                :
                    tracked_(true),
                    inPatchId_(inPatchId),
                    storePositions_(false),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    currentPosition_(initialPosition),
                    distanceTravelled_(vector::zero),
                    msd_(vector::zero)//,
                    //parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const label& inPatchId,
                    const scalar& initialTime,
                    const vector& initialPosition
                )
                :
                    tracked_(tracked),
                    inPatchId_(inPatchId),
                    storePositions_(false),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    currentPosition_(initialPosition),
                    distanceTravelled_(vector::zero),
                    msd_(vector::zero)//,
                    //parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const label& inPatchId,
                    const bool& storePositions,
                    const scalar& initialTime,
                    const vector& initialPosition
                )
                :
                    tracked_(tracked),
                    inPatchId_(inPatchId),
                    storePositions_(storePositions),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    currentPosition_(initialPosition),
                    distanceTravelled_(vector::zero),
                    msd_(vector::zero)//,
                    //parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const label& inPatchId,
                    const scalar& initialTime,
                    const vector& initialPosition,
                    const vector& distanceTravelled
                )
                :
                    tracked_(tracked),
                    inPatchId_(inPatchId),
                    storePositions_(false),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    currentPosition_(initialPosition),
                    distanceTravelled_(distanceTravelled),
                    msd_(vector::zero)//,
                    //parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const label& inPatchId,
                    const bool& storePositions,
                    const scalar& initialTime,
                    const vector& initialPosition,
                    const vector& distanceTravelled//,
                    //const DynamicList<FixedList<scalar, 4>>& parcelTrajectory
                )
                :
                    tracked_(tracked),
                    inPatchId_(inPatchId),
                    storePositions_(storePositions),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    currentPosition_(initialPosition),
                    distanceTravelled_(distanceTravelled),
                    msd_(vector::zero)//,
                    //parcelTrajectory_(parcelTrajectory)
                {}
                
                //- Copy constructor
                TrackedParcel
                (
                    const TrackedParcel& tP
                )
                :
                    tracked_(tP.tracked_),
                    inPatchId_(tP.inPatchId_),
                    storePositions_(tP.storePositions_),
                    t0_(tP.t0_),
                    initialPosition_(tP.initialPosition_),
                    currentPosition_(tP.currentPosition_),
                    distanceTravelled_(tP.distanceTravelled_),
                    msd_(tP.msd_)//,
                    //parcelTrajectory_(tP.parcelTrajectory_)
                {}
                
                //- Construct from Istream
                TrackedParcel(Istream& is)
                :
                    tracked_(false),
                    inPatchId_(-1),
                    storePositions_(false),
                    t0_(0.0),
                    initialPosition_(vector::zero),
                    currentPosition_(vector::zero),
                    distanceTravelled_(vector::zero),
                    msd_(vector::zero)//,
                    //parcelTrajectory_()
                {
                    Info << "Construct from Istream obj" << endl;
                    if (is.format() == IOstream::ASCII)
                    {
                        tracked_ = readBool(is);
                        inPatchId_ = readLabel(is);
                        storePositions_ = readBool(is);
                        t0_ = readScalar(is);
                        is >> initialPosition_
                           >> currentPosition_ 
                           >> distanceTravelled_ 
                           >> msd_;
                        //is >> parcelTrajectory_;
                    }
                    else
                    {
                        is.read
                        (
                            reinterpret_cast<char*>(&initialPosition_),
                            sizeof(tracked_)
                            + sizeof(inPatchId_)
                            + sizeof(storePositions_)
                            + sizeof(t0_)
                            + sizeof(initialPosition_)
                            //+ sizeof(currentPosition_) // careful
                            //+ sizeof(distanceTravelled_) // careful
                            //+ sizeof(msd_) // careful
                        );
                        
                        is >> currentPosition_ 
                           >> distanceTravelled_ 
                           >> msd_;
                        
                        //is>> parcelTrajectory_;
                    }
                    
                    // Check state of Istream
                    is.check("TrackedParcel::TrackedParcel(Istream&)");
                }
                
            
            // Destructor
                
                ~TrackedParcel(){}
                
            
            // Member functions
            
                //- Set the initial parcel info
                void setInitialParcelInfo
                (
                    const label& inPatchId,
                    const scalar& t0,
                    const vector& iniPos
                )
                {
                    inPatchId_ = inPatchId;
                    t0_ = t0;
                    initialPosition_ = iniPos;
                    currentPosition_ = initialPosition_;
                }
                
                //- Set the initial position of the parcel
                void setInitialParcelPosition
                (
                    const scalar& t0 = 0.0,
                    const vector& iniPos = vector::zero
                )
                {
                    t0_ = t0;
                    initialPosition_ = iniPos;
                    currentPosition_ = initialPosition_;
                }
                
                //- Update the parcel's trajectory
                void updateParcelTrajectory(const scalar t, const vector& pos)
                {
                    FixedList<scalar, 4> sF(0.0);
                    sF[0] = t;
                    sF[1] = pos[0];
                    sF[2] = pos[1];
                    sF[3] = pos[2];
                    //parcelTrajectory_.append(sF);
                }
            
            
            // Access
            
                //- Return const access to tracked
                inline const bool& tracked() const
                {
                    return tracked_;
                }
                
                //- Return const access to inPatchId
                inline label inPatchId() const
                {
                    return inPatchId_;
                };
                
                //- Return const access to storePositions
                inline const bool& storePositions() const
                {
                    return storePositions_;
                }
                
                //- Return const access to the parcel's initial time
                inline const scalar& initialTime() const
                {
                    return t0_;
                }
                
                //- Return const access to the parcel's initial position
                inline const vector& initialPosition() const
                {
                    return initialPosition_;
                }
                
                //- Return const access to the parcel's current position
                inline const vector& currentPosition() const
                {
                    return currentPosition_;
                }
                
                //- Return const access to the distance travelled by the parcel
                inline scalar distanceTravelled() const
                {
                    return mag(distanceTravelled_);
                }
                
                //- Return const access to the parcel's mean squared 
                //  displacement
                inline scalar meanSquaredDisplacement() const
                {
                    return msd_[0] + msd_[1] + msd_[2];
                }
                
                //- Return const access to the distance travelled by the parcel
                //  in vector form
                inline const vector& distanceTravelledVector() const
                {
                    return distanceTravelled_;
                }
                
                //- Return const access to the parcel's mean squared 
                //  displacement in vector form
                inline const vector& meanSquaredDisplacementVector() const
                {
                    return msd_;
                }
                
                //- Return const access to the full parcel's trajectory
                /*inline const DynamicList<FixedList<scalar, 4>>& 
                parcelTrajectory() const
                {
                    return parcelTrajectory_;
                }*/
            
            
            // Edit
            
                //- Return access to tracked
                inline bool& tracked()
                {
                    return tracked_;
                }
                
                //- Return access to inPatchId
                inline label& inPatchId()
                {
                    return inPatchId_;
                };
                
                //- Return access to storePositions
                inline bool& storePositions()
                {
                    return storePositions_;
                }
                
                //- Return access to the parcel's initial time
                inline scalar& initialTime()
                {
                    return t0_;
                }
                
                //- Return access to the parcel's initial position
                inline vector& initialPosition()
                {
                    return initialPosition_;
                }
                
                //- Return access to the parcel's current position
                inline vector& currentPosition()
                {
                    return currentPosition_;
                }
                
                //- Update the current position of the parcel
                void updateCurrentPosition(const vector& currentPos)
                {
                    currentPosition_ = currentPos;
                }
                
                //- Return access to the distance travelled by the parcel in
                //  vector form
                inline vector& distanceTravelledVector()
                {
                    return distanceTravelled_;
                }
                
                //- Update the distance travelled by the parcel
                void updateDistanceTravelled(const vector& newPos)
                {
                    distanceTravelled_ += newPos - currentPosition_;
                    
                    updateCurrentPosition(newPos);
                }
                
                //- Update the total distance travelled by the parcel
                void updateTotalDistanceTravelled(const vector& newPos)
                {
                    distanceTravelled_ += vector
                        (
                            mag(newPos[0] - currentPosition_[0]), 
                            mag(newPos[1] - currentPosition_[1]), 
                            mag(newPos[2] - currentPosition_[2])
                        );
                        
                    updateCurrentPosition(newPos);
                }
                
                //- Update the distance travelled by the parcel
                void performSpecularReflectionOnDistanceTravelled
                (
                    const vector& nw
                )
                {
                    forAll(nw, dim)
                    {
                        if(abs(nw[dim]) > SMALL)
                        {
                            distanceTravelled_[dim] *= -1;
                        }
                    }
                }
                
                //- Return access to the parcel's mean squared displacement in
                //  vector form
                inline vector& meanSquaredDisplacementVector()
                {
                    return msd_;
                }
                
                //- Update the parcel's mean squared displacement
                void updateMeanSquaredDisplacement()
                {
                    msd_ = vector
                        (
                            sqr(distanceTravelled_[0]), 
                            sqr(distanceTravelled_[1]), 
                            sqr(distanceTravelled_[2])
                        ); 
                }
                
                //- Return access to the parcel's full trajectory
                /*inline DynamicList<FixedList<scalar, 4>>& parcelTrajectory()
                {
                    return parcelTrajectory_;
                }*/
                
                
            // Operators
                
                //- Assignment operator
                TrackedParcel& operator=(const TrackedParcel& tP)
                {
                     tracked_ = tP.tracked_;
                     inPatchId_ = tP.inPatchId_;
                     storePositions_ = tP.storePositions_;
                     t0_ = tP.t0_;
                     initialPosition_ = tP.initialPosition_;
                     currentPosition_ = tP.currentPosition_;
                     distanceTravelled_ = tP.distanceTravelled_;
                     msd_ = tP.msd_;
                     //parcelTrajectory_ = tP.parcelTrajectory_;
                     
                     return *this;
                }
                
                //- Istream operator
                friend Istream& operator>>         
                (
                    Istream& is,
                    TrackedParcel& tP
                )
                {
                    if (is.format() == IOstream::ASCII)
                    {
                        tP.tracked_ = readBool(is);
                        tP.inPatchId_ = readLabel(is);
                        tP.storePositions_ = readBool(is);
                        tP.t0_ = readScalar(is);
                        is >> tP.initialPosition_ 
                           >> tP.currentPosition_
                           >> tP.distanceTravelled_ 
                           >> tP.msd_;
                        //is >> tP.parcelTrajectory_;
                    }
                    else
                    {
                        is.read
                        (
                            reinterpret_cast<char*>(&tP.initialPosition_),
                            sizeof(tP.tracked())
                            + sizeof(tP.inPatchId())
                            + sizeof(tP.storePositions())
                            + sizeof(tP.initialTime())
                            + sizeof(tP.initialPosition())
                            //+ sizeof(tP.currentPosition()) // careful
                            //+ sizeof(tP.distanceTravelledVector()) // careful
                            //+ sizeof(tP.meanSquaredDisplacementVector()) // careful
                        );
                        
                        is >> tP.currentPosition_
                           >> tP.distanceTravelled_
                           >> tP.msd_;
                        
                        //is >> tP.parcelTrajectory_;
                    }
                    
                    // Check state of Istream
                    is.check
                    (
                        "Foam::Istream& Foam::operator>>"
                        "(Foam::Istream&, Foam::dsmcParcel::TrackedParcel&)"
                    );
                    
                    return is;
                }
                
                //- Ostream operator
                friend Ostream& operator<<         
                (
                    Ostream& os,
                    const TrackedParcel& tP
                )
                {
                    if (os.format() == IOstream::ASCII)
                    {
                        os << tP.tracked_ << token::SPACE
                           << tP.inPatchId_ << token::SPACE
                           << tP.storePositions_ << token::SPACE
                           << tP.t0_ << token::SPACE
                           << tP.initialPosition_ << token::SPACE
                           << tP.currentPosition_ << token::SPACE
                           << tP.distanceTravelled_ << token::SPACE
                           << tP.msd_ //<< token::SPACE
                           /*<< tP.parcelTrajectory_*/;
                    }
                    else
                    {
                        os.write
                        (
                            reinterpret_cast<const char*>
                                (&tP.initialPosition_),
                            sizeof(tP.tracked())
                            + sizeof(tP.inPatchId())
                            + sizeof(tP.storePositions())
                            + sizeof(tP.initialTime())
                            + sizeof(tP.initialPosition())
                            //+ sizeof(tP.currentPosition()) // careful
                            //+ sizeof(tP.distanceTravelledVector()) // careful
                            //+ sizeof(tP.meanSquaredDisplacementVector()) // careful
                        );
                        
                        os << tP.currentPosition_
                           << tP.distanceTravelled_
                           << tP.msd_;
                           
                        //os << tP.parcelTrajectory_;
                    }
                       
                    // Check state of Ostream
                    os.check
                    (
                        "Foam::Ostream& Foam::operator<<"
                        "(Foam::Ostream&, const Foam::dsmcParcel::TrackedParcel&)"
                    );
                    
                    return os;
                }
    };
    
    
    //- Class used to pass kinematic tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<dsmcCloud>
    {
        // Private data



    public:

        // Constructors

            //- Construct from components
            trackingData(dsmcCloud& cloud)
            :
                particle::TrackingData<dsmcCloud>(cloud)
            {}

    };


protected:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;
            
            //- Radial weight of parcel (for axisymmetric)
            scalar RWF_;

            //- Rotational energy of the Parcel [J]
            scalar ERot_;
            
            //- Electronic energy level of the Parcel
            label ELevel_;

            //- Parcel type id
            label typeId_;
            
            //- Patch iD of newly inserted parcel (default: -1)
            label newParcel_;
            
            //- Pointer to TrackedParcel object containing tracking info if the
            //  parcel is being selected -> TO BE REINTRODUCED
            //TrackedParcel* tracked_;
            //- TrackedParcel object containing tracking info if the
            //  parcel is being selected -> TEMPORARY
            TrackedParcel tracked_;
            
            //- see W.F.N. Santos, "Physical and computational aspects
            //  of shock waves over power-law leading edges",
            //  Physics of Fluids, 20, 016101, 2008.
            label classification_;
            
            //- Pointer to StuckParcel object containing stuck info if the
            //  parcel has been adsorbed
            StuckParcel* stuck_;
            
            //- Vibrational energy level for each vibrational energy mode of
            //  the Parcel
            labelList vibLevel_;

 
public:

    //- Runtime type information
//     TypeName("dsmcParcel");

    friend class Cloud<dsmcParcel>;


    // Constructors

        //- Construct from components
        inline dsmcParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& U,
            const scalar RWF,
            const scalar ERot,
            const label ELevel,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const label newParcel,
            const label classification,
            const labelList& vibLevel
        );
        
        //- Construct from Istream
        /*dsmcParcel
          (
             const Cloud<dsmcParcel>& c,
             Istream& is,
             bool readFields = true
          );*/
        
        dsmcParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<dsmcParcel> clone() const
        {
            return autoPtr<dsmcParcel>(new dsmcParcel(*this));
        }
        
        //- Copy constructor
        dsmcParcel
        (
            const dsmcParcel& dP
        )
        :
            particle(dP),
            U_(dP.U_),
            RWF_(dP.RWF_),
            ERot_(dP.ERot_),
            ELevel_(dP.ELevel_),
            typeId_(dP.typeId_),
            newParcel_(dP.newParcel_),
            //tracked_(nullptr), // TODO TO BE REINTRODUCED
            tracked_(dP.tracked_),
            classification_(dP.classification_),
            stuck_(nullptr),
            vibLevel_(dP.vibLevel_)
        {
            /*if(dP.tracked_)
            {
                setTracked();
                *tracked_ = *(dP.tracked_);
            }*/ // TODO TO BE REINTRODUCED
            
            if(dP.stuck_)
            {
                setStuck();
                *stuck_ = *(dP.stuck_);
            }
        }
        
        
    // Destructor
                
        ~dsmcParcel()
        {
            //deleteTracked(); // TODO TO BE REINTRODUCED
            
            deleteStuck();
        }



        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<dsmcParcel> operator()(Istream& is) const
            {
                return autoPtr<dsmcParcel>(new dsmcParcel(mesh_, is, true));
            }
        };

        
    // Member Functions

        // Access

            //- Return const access to type id
            inline const label& typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;
            
            //- Return const access to radial weight
            inline const scalar& RWF() const;

            //- Return const access to rotational energy
            inline const scalar& ERot() const;
            
            //- Return const access to electronic energy
            inline const label& ELevel() const;
            
            //- Return const access to new parcel
            inline const label& newParcel() const;
            
            //- Return const access to the TrackedParcel object pointed by 
            //  tracked_
            inline const TrackedParcel& tracked() const;
            
            //- Return the state of the tracked_ pointer
            //  True if active, false otherwise
            inline bool isTracked() const;
            
            //- Return const access to parcel classification 
            //  see Santos W.F., Physical and computational 
            //  aspects of shock waves over power-law leading edges,
            //  Physics of Fluids, 20, 016101, 2008
            inline const label& classification() const;
            
            //- Return const access to the StuckParcel object pointed by 
            //  stuck_
            inline const StuckParcel& stuck() const;
            
            //- Return the state of the stuck_ pointer
            //  True if active, false otherwise
            inline bool isStuck() const;
            
            //- Return the state of the stuck_ pointer
            //  False if active, true otherwise
            inline bool isFree() const;
            
            //- Return const access to vibrational energy level for each vib
            //  energy mode
            inline const labelList& vibLevel() const;


        // Edit

            //- Return non-const access to typeId
            inline label& typeId();

            //- Return access to velocity
            inline vector& U();
            
            //- Return non-const access to radial weight
            inline scalar& RWF();

            //- Return access to rotational energy
            inline scalar& ERot();
            
            //- Return access to electronic energy
            inline label& ELevel();
            
            //- Return access to new parcel
            inline label& newParcel();
            
            //- Return access to the TrackedParcel object
            inline TrackedParcel& tracked();
            
            //- Construct the TrackedParcel object pointed by tracked_
            inline void setTracked();
            
            inline void setTracked
            (
                const bool& tracked,
                const label& inPatchId,
                const scalar& initialTime,
                const vector& initialPosition
            );
            
            inline void setTracked
            (
                const bool& tracked,
                const label& inPatchId,
                const bool& storePositions,
                const scalar& initialTime,
                const vector& initialPosition
            );
            
            inline void setTracked
            (
                const bool& tracked,
                const label& inPatchId,
                const scalar& initialTime,
                const vector& initialPosition,
                const vector& distanceTravelled
            );
            
            inline void setTracked
            (
                const bool& tracked,
                const label& inPatchId,
                const bool& storePositions,
                const scalar& initialTime,
                const vector& initialPosition,
                const vector& distanceTravelled/*,
                const DynamicList<FixedList<scalar, 4>>& parcelTrajectory*/
            );
            
            //- Delete the TrackedParcel pointer and set it to NULL
            inline void deleteTracked();
            
            //- Return access to parcel classification 
            //  see Santos W.F., Physical and computational 
            //  aspects of shock waves over power-law leading edges,
            //  Physics of Fluids, 20, 016101, 2008
            inline label& classification();
            
            //- Return access to the StuckParcel object
            inline StuckParcel& stuck();
            
            //- Construct the StuckParcel object pointed by stuck_
            inline void setStuck();
            
            inline void setStuck
            (
                const scalarField& wallTemperature,
                const vectorField& wallVectors
            );
            
            //- Delete the StuckParcel pointer and set it to NULL
            inline void deleteStuck();
            
            //- Return access to vibrational energy
            inline labelList& vibLevel();


        // Main calculation loop

            // Tracking

                //- Move the parcel
                bool move(trackingData& td, const scalar trackTime);
                

        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            //  template<class trackData>
            //  Patch interactions

            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackingData& td,
                const tetIndices&
            );
            
            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );
            
            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties(const vector& separation);


        // I-O

            static void readFields(Cloud<dsmcParcel>& c);

            static void writeFields(const Cloud<dsmcParcel>& c);


    // Operators
    
        //- Assignment operator
        dsmcParcel& operator=(const dsmcParcel& dP)
        {
            particle::operator=(dP);
            U_ = dP.U_;
            RWF_ = dP.RWF_;
            ERot_ = dP.ERot_;
            ELevel_ = dP.ELevel_;
            typeId_ = dP.typeId_;
            newParcel_ = dP.newParcel_;
            
            tracked_ = dP.tracked_;
            /*if (dP.isTracked())
            {
                setTracked();
                *tracked_ = *(dP.tracked_);
            }
            else
            {
                tracked_ = nullptr;
            }*/ // TODO TO BE REINTRODUCED
            
            classification_ = dP.classification_;
            
            if (dP.isStuck())
            {
                setStuck();
                *stuck_ = *(dP.stuck_);
            }
            else
            {
                stuck_ = nullptr;
            }
            
            vibLevel_ = dP.vibLevel_;

            return *this;
        }
        
        
    // Friend Operators
        
        //- Ostream operator
        friend Ostream& operator<<         
        (
            Ostream&,
            const dsmcParcel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
