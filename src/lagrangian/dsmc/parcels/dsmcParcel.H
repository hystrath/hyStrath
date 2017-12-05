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

            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Rotational degrees of freedom
            scalar rotationalDegreesOfFreedom_;
        
            //- Vibrational degrees of freedom
            scalar vibrationalDegreesOfFreedom_;

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
            
            //- Reference temperature for reference vibrational relaxation number
            scalarList TrefZv_;
            
            //- Characteristic dissociation quantum level
            labelList charDissQuantumLevel_;
            
            //- Ionisation temperature
            scalar ionisationTemperature_;
            
            //- Number of electronic energy/degeneracy levels
            label numberOfElectronicLevels_;
            
            //- Charge state (0 for neutrals, 1 for ions, -1 for electrons)
            label charge_;

            //- Scalar list of electronic energy levels [J]
            scalarList electronicEnergyList_;     
        
            //- Integer list of degeneracies at each electronic level
            labelList degeneracyList_;

    public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member functions

            //- Return const access to the particle density
            inline scalar mass() const;

            //- Return const access to the minimum particle mass
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the rotationalDegreesOfFreedom
            inline scalar rotationalDegreesOfFreedom() const;
            
            //- Return the vibrationallDegreesOfFreedom
            inline scalar vibrationalDegreesOfFreedom() const;
            
            //- Return the viscosity index
            inline scalar omega() const;
            
            //- Return the VSS scattering parameter
            inline scalar alpha() const;
    
            //- Return the characteristic vibrational temperature
            inline scalarList thetaV() const;    
            
            //- Return the characteristic dissociation temperature
            inline scalarList thetaD() const;
    
            //- Return the characteristic vibrational relaxation number
            inline scalarList Zref() const;
                      
            //- Return the reference temperature for the characteristic vibrational relaxation number
            inline scalarList TrefZv() const;
            
            //- Return the characteristic dissociation level
            inline labelList charDissQuantumLevel() const;
            
            //- Return the ionisation temperature
            inline scalar ionisationTemperature() const;
            
            //- Return the number of electronic energy/degeneracy levels
            inline label numberOfElectronicLevels() const;
            
            //- Return the particle's charge state
            inline label charge() const;

            //- Return the number of electronic energy/degeneracy levels
            inline scalarList electronicEnergyList() const;
            
            //- Return the number of electronic energy/degeneracy levels
            inline labelList degeneracyList() const;
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
              
              //- Effective diffusivity for all tracked parcels
              static scalar D_EFF;
        
        protected:
        
            // Protected data
            
              //- Is the parcel being tracked
              bool tracked_;
              
              //- Record the parcel's trajectory
              bool storePositions_;
              
              //- Time at which the parcel was selected as a tracer
              scalar t0_;
              
              //- Parcel's initial position
              vector initialPosition_;
              
              //- Mean square displacement
              vector msd_;
              
              //- Effective diffusivity
              scalar deff_;
              
              //- A data structure holding the position of the tracked parcel
              //  versus time
              DynamicList<FixedList<scalar, 4>> parcelTrajectory_;
            
            
            // Protected member functions
            
              //- Calculation of the parcel's mean square displacement
              void msd(const vector& currentPos)
              {
                  const vector& diffPos = currentPos - initialPosition_;
                    
                  msd_ = vector
                      (
                          sqr(diffPos.component(0)), 
                          sqr(diffPos.component(1)), 
                          sqr(diffPos.component(2))
                      );
              }
              
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
                    storePositions_(false),
                    t0_(0.0),
                    initialPosition_(vector::zero),
                    msd_(vector::zero),
                    deff_(0.0),
                    parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const scalar& initialTime,
                    const vector& initialPosition
                )
                :
                    tracked_(true),
                    storePositions_(false),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    msd_(vector::zero),
                    deff_(0.0),
                    parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const scalar& initialTime,
                    const vector& initialPosition
                )
                :
                    tracked_(tracked),
                    storePositions_(false),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    msd_(vector::zero),
                    deff_(0.0),
                    parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const bool& storePositions,
                    const scalar& initialTime,
                    const vector& initialPosition
                )
                :
                    tracked_(tracked),
                    storePositions_(storePositions),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    msd_(vector::zero),
                    deff_(0.0),
                    parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const scalar& initialTime,
                    const vector& initialPosition,
                    const vector& msd
                )
                :
                    tracked_(tracked),
                    storePositions_(false),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    msd_(msd),
                    deff_(0.0),
                    parcelTrajectory_()
                {}
                
                TrackedParcel
                (
                    const bool& tracked,
                    const bool& storePositions,
                    const scalar& initialTime,
                    const vector& initialPosition,
                    const vector& msd,
                    const DynamicList<FixedList<scalar, 4>>& parcelTrajectory
                )
                :
                    tracked_(tracked),
                    storePositions_(storePositions),
                    t0_(initialTime),
                    initialPosition_(initialPosition),
                    msd_(msd),
                    deff_(0.0),
                    parcelTrajectory_(parcelTrajectory)
                {}
                
                //- Copy constructor
                TrackedParcel
                (
                    const TrackedParcel& tP
                )
                :
                    tracked_(tP.tracked_),
                    storePositions_(tP.storePositions_),
                    t0_(tP.t0_),
                    initialPosition_(tP.initialPosition_),
                    msd_(tP.msd_),
                    deff_(tP.deff_),
                    parcelTrajectory_(tP.parcelTrajectory_)
                {}
                
                //- Construct from Istream
                TrackedParcel(Istream& is)
                :
                    tracked_(false),
                    storePositions_(false),
                    t0_(0.0),
                    initialPosition_(vector::zero),
                    msd_(vector::zero),
                    deff_(0.0),
                    parcelTrajectory_()
                {
                    if (is.format() == IOstream::ASCII)
                    {
                        tracked_ = readBool(is);
                        storePositions_ = readBool(is);
                        t0_ = readScalar(is);
                        is >> initialPosition_ >> msd_;
                        deff_ = readScalar(is);
                        is >> parcelTrajectory_;
                    }
                    else
                    {
                        is.read
                        (
                            reinterpret_cast<char*>(&initialPosition_),
                            sizeof(tracked_)
                            + sizeof(storePositions_)
                            + sizeof(t0_)
                            + sizeof(initialPosition_)
                            + sizeof(deff_)
                        );
                        
                        is >> parcelTrajectory_;
                    }
                    
                    // Check state of Istream
                    is.check("TrackedParcel::TrackedParcel(Istream&)");
                }
                
            
            // Destructor
                
                ~TrackedParcel()
                {
                    if(msd_[0] >= 0.99)
                    {
                        nDELETED++;
                    
                        D_EFF += deff_;
                        
                        Info << "deff porous medium: " << D_EFF/nDELETED << endl;
                    }
                }

                
            
            // Member functions
            
                //- Set the initial parcel info
                void setInitialParcelInfo
                (
                    const bool& tracked = true,
                    const scalar& t0 = 0.0,
                    const vector& iniPos = vector::zero,
                    const vector& msd = vector::zero
                )
                {
                    tracked_ = tracked;
                    t0_ = t0;
                    initialPosition_ = iniPos;
                    msd_ = msd;
                }
                
                //- Set the initial position of the parcel
                void setInitialParcelPosition
                (
                    const vector& iniPos = vector::zero
                )
                {
                    initialPosition_ = iniPos;
                }
                
                //- Update the parcel's trajectory
                void updateParcelTrajectory(const scalar t, const vector& pos)
                {
                    FixedList<scalar, 4> sF(0.0);
                    sF[0] = t;
                    sF[1] = pos[0];
                    sF[2] = pos[1];
                    sF[3] = pos[2];
                    parcelTrajectory_.append(sF);
                }
            
            
            // Access
            
                //- Return const access to tracked
                inline const bool& tracked() const
                {
                    return tracked_;
                }
                
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
                
                //- Return the parcel's mean square displacement
                inline scalar meanSquareDisplacement() const
                {
                    return msd_[0] + msd_[1] + msd_[2];
                }
                
                //- Return the parcel's mean square displacement
                //  given the parcel's current location
                scalar meanSquareDisplacement(const vector& currentPosition)
                {
                    msd(currentPosition);
                    
                    return meanSquareDisplacement();
                }
                
                //- Return const access to the parcel's mean square 
                //  displacement in vector form
                inline const vector& meanSquareDisplacementVector() const
                {
                    return msd_;
                }
                
                //- Return const access to the parcel's effective diffusivity
                inline const scalar& deff() const
                {
                    return deff_;
                }
                
                //- Return const access to the full parcel's trajectory
                inline const DynamicList<FixedList<scalar, 4>>& 
                parcelTrajectory() const
                {
                    return parcelTrajectory_;
                }
            
            
            // Edit
            
                //- Return access to tracked
                inline bool& tracked()
                {
                    return tracked_;
                }
                
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
                
                //- Return access to the parcel's mean square displacement in
                //  vector form
                inline vector& meanSquareDisplacementVector()
                {
                    return msd_;
                }
                
                //- Return access to the parcel's effective diffusivity
                inline scalar& deff()
                {
                    return deff_;
                }
                
                //- Return access to the parcel's full trajectory
                inline DynamicList<FixedList<scalar, 4>>& parcelTrajectory()
                {
                    return parcelTrajectory_;
                }
                
                
            // Operators
                
                //- Assignment operator
                TrackedParcel& operator=(const TrackedParcel& tP)
                {
                     tracked_ = tP.tracked_;
                     storePositions_ = tP.storePositions_;
                     t0_ = tP.t0_;
                     initialPosition_ = tP.initialPosition_;
                     msd_ = tP.msd_;
                     deff_ = tP.deff_;
                     parcelTrajectory_ = tP.parcelTrajectory_;
                     
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
                        tP.storePositions_ = readBool(is);
                        tP.t0_ = readScalar(is);
                        is >> tP.initialPosition_;
                        tP.deff_ = readScalar(is);
                        is >> tP.parcelTrajectory_;
                    }
                    else
                    {
                        is.read
                        (
                            reinterpret_cast<char*>(&tP.initialPosition_),
                            sizeof(tP.tracked())
                            + sizeof(tP.storePositions())
                            + sizeof(tP.initialTime())
                            + sizeof(tP.initialPosition())
                            + sizeof(tP.deff_)
                        );
                        
                        is >> tP.parcelTrajectory_;
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
                           << tP.storePositions_ << token::SPACE
                           << tP.t0_ << token::SPACE
                           << tP.initialPosition_ << token::SPACE
                           << tP.deff_ << token::SPACE
                           << tP.parcelTrajectory_;
                    }
                    else
                    {
                        os.write
                        (
                            reinterpret_cast<const char*>
                                (&tP.initialPosition_),
                            sizeof(tP.tracked())
                            + sizeof(tP.storePositions())
                            + sizeof(tP.initialTime())
                            + sizeof(tP.initialPosition())
                             + sizeof(tP.deff_)
                        );
                        
                        os << tP.parcelTrajectory_;
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


private:

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
            
            //- newly inserted parcel
            label newParcel_;
            
            //- Pointer to TrackedParcel object containing tracking info if the
            //  parcel is being selected
            TrackedParcel* tracked_;
            
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
            tracked_(nullptr),
            classification_(dP.classification_),
            stuck_(nullptr),
            vibLevel_(dP.vibLevel_)
        {
            if(dP.tracked_)
            {
                setTracked();
                *tracked_ = *(dP.tracked_);
            }
            
            if(dP.stuck_)
            {
                setStuck();
                *stuck_ = *(dP.stuck_);
            }
        }
        
        
    // Destructor
                
        ~dsmcParcel()
        {
            delete tracked_;
            tracked_ = nullptr;
            
            delete stuck_;
            stuck_ = nullptr;
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
                const scalar& initialTime,
                const vector& initialPosition
            );
            
            inline void setTracked
            (
                const bool& tracked,
                const bool& storePositions,
                const scalar& initialTime,
                const vector& initialPosition
            );
            
            inline void setTracked
            (
                const bool& tracked,
                const scalar& initialTime,
                const vector& initialPosition,
                const vector& meanSquareDisplacement
            );
            
            //- Construct the TrackedParcel object pointed by tracked_
            inline void setTracked
            (
                const bool& tracked,
                const bool& storePositions,
                const scalar& initialTime,
                const vector& initialPosition,
                const vector& meanSquareDisplacement,
                const DynamicList<FixedList<scalar, 4>>& parcelTrajectory
            );
            
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
            
            if (dP.isTracked())
            {
                setTracked();
                *tracked_ = *(dP.tracked_);
            }
            else
            {
                tracked_ = nullptr;
            }
            
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
