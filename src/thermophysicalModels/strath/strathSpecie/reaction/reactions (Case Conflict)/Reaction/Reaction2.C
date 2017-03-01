/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "Reaction2.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

template<class Reaction2Thermo>
Foam::label Foam::Reaction2<Reaction2Thermo>::nUnNamedReactions = 0;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Reaction2Thermo>
void Foam::Reaction2<Reaction2Thermo>::reactionStrLeft
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < lhs_.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(lhs_[i].stoichCoeff - 1) > SMALL)
        {
            reaction << lhs_[i].stoichCoeff;
        }
        reaction << species_[lhs_[i].index];
        if (mag(lhs_[i].exponent - lhs_[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << lhs_[i].exponent;
        }
    }
}


template<class Reaction2Thermo>
void Foam::Reaction2<Reaction2Thermo>::reactionStrRight
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < rhs_.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(rhs_[i].stoichCoeff - 1) > SMALL)
        {
            reaction << rhs_[i].stoichCoeff;
        }
        reaction << species_[rhs_[i].index];
        if (mag(rhs_[i].exponent - rhs_[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << rhs_[i].exponent;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Reaction2Thermo>
Foam::label Foam::Reaction2<Reaction2Thermo>::getNewReactionID()
{
    return nUnNamedReactions++;
}


template<class Reaction2Thermo>
Foam::string Foam::Reaction2<Reaction2Thermo>::reactionStr
(
    OStringStream& reaction
) const
{
    reactionStrLeft(reaction);
    reaction << " = ";
    reactionStrRight(reaction);
    return reaction.str();
}


template<class Reaction2Thermo>
void Foam::Reaction2<Reaction2Thermo>::setThermo
(
    const HashPtrTable<Reaction2Thermo>& thermoDatabase
)
{
    if (rhs_.size() > 0)
    {
        Reaction2Thermo::operator=
        (
            rhs_[0].stoichCoeff*(*thermoDatabase[species_[rhs_[0].index]])
        );

        for (label i=1; i<rhs_.size(); ++i)
        {
            this->operator+=
            (
                rhs_[i].stoichCoeff*(*thermoDatabase[species_[rhs_[i].index]])
            );
        }
    }

    forAll(lhs_, i)
    {
        this->operator-=
        (
            lhs_[i].stoichCoeff*(*thermoDatabase[species_[lhs_[i].index]])
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Reaction2Thermo>
Foam::Reaction2<Reaction2Thermo>::Reaction2
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const HashPtrTable<Reaction2Thermo>& thermoDatabase
)
:
    Reaction2Thermo(*thermoDatabase[species[0]]),
    name_("un-named-reaction-" + Foam::name(getNewReactionID())),
    species_(species),
    lhs_(lhs),
    rhs_(rhs)
{
    // NEW VINCENT - UNUSED CONSTRUCTOR
    word controllingTemperature = "dissociation";
    if (controllingTemperature == "dissociation")
    {
        controlT_ = dissociation;
    }
    Info << "This constructor should not be used - VINCENT" << endl;
    // END NEW VINCENT
    setThermo(thermoDatabase);
}


template<class Reaction2Thermo>
Foam::Reaction2<Reaction2Thermo>::Reaction2
(
    const Reaction2<Reaction2Thermo>& r,
    const speciesTable& species
)
:
    Reaction2Thermo(r),
    name_(r.name() + "Copy"),
    species_(species),
    lhs_(r.lhs_),
    rhs_(r.rhs_),
    controlT_(r.controlT_) // NEW VINCENT
{}


template<class Reaction2Thermo>
Foam::Reaction2<Reaction2Thermo>::specieCoeffs::specieCoeffs
(
    const speciesTable& species,
    Istream& is
)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1.0;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        size_t i = specieName.find('^');

        if (i != word::npos)
        {
            string exponentStr = specieName
            (
                i + 1,
                specieName.size() - i - 1
            );
            exponent = atof(exponentStr.c_str());
            specieName = specieName(0, i);
        }

        if (species.contains(specieName))
        {
            index = species[specieName];
        }
        else
        {
            index = -1;
        }
    }
    else
    {
        FatalIOErrorIn("Reaction2<Reaction2Thermo>::lrhs(Istream& is)", is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}


template<class Reaction2Thermo>
void Foam::Reaction2<Reaction2Thermo>::setLRhs
(
    Istream& is,
    const speciesTable& species,
    List<specieCoeffs>& lhs,
    List<specieCoeffs>& rhs
)
{
    DynamicList<specieCoeffs> dlrhs;

    while (is.good())
    {
        dlrhs.append(specieCoeffs(species, is));

        if (dlrhs.last().index != -1)
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    lhs = dlrhs.shrink();
                    dlrhs.clear();
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                rhs = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            dlrhs.remove();
            if (is.good())
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
            }
            else
            {
                if (!dlrhs.empty())
                {
                    rhs = dlrhs.shrink();
                }
                return;
            }
        }
    }

    FatalIOErrorIn("Reaction2<Reaction2Thermo>::setLRhs(Istream& is)", is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


template<class Reaction2Thermo>
Foam::Reaction2<Reaction2Thermo>::Reaction2
(
    const speciesTable& species,
    const HashPtrTable<Reaction2Thermo>& thermoDatabase,
    Istream& is
)
:
    Reaction2Thermo(*thermoDatabase[species[0]]),
    name_("un-named-reaction" + Foam::name(getNewReactionID())),
    species_(species)
{
    Info << "This constructor should not be used - VINCENT" << endl;
    // NEW VINCENT
    word controllingTemperature = "dissociation";
    if (controllingTemperature == "dissociation")
    {
        controlT_ = dissociation;
    }
    // END NEW VINCENT
    
    setLRhs(is, species, lhs_, rhs_);
    setThermo(thermoDatabase);
}


template<class Reaction2Thermo>
Foam::Reaction2<Reaction2Thermo>::Reaction2
(
    const speciesTable& species,
    const HashPtrTable<Reaction2Thermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction2Thermo(*thermoDatabase[species[0]]),
    name_(dict.dictName()),
    species_(species)
{
    // This constructor is used in the very beginnning of a simulation
    // BRAND NEW VINCENT ******************************************************
    word controllingTemperature = dict.lookupOrDefault<word>("controlT", "transrotational");
    if (controllingTemperature == "chargeExchange")
    {
        controlT_ = chargeExchange;
    }
    else if (controllingTemperature == "dissociation")
    {
        controlT_ = dissociation;
    }
    else if (controllingTemperature == "exchange")
    {
        controlT_ = exchange;
    }
    else if (controllingTemperature == "impactDissociation")
    {
        controlT_ = impactDissociation;
    }
    else if (controllingTemperature == "impactIonisation")
    {
        controlT_ = impactIonisation;
    }
    else if (controllingTemperature == "associativeIonisation")
    {
        controlT_ = associativeIonisation;
    }
    else if (controllingTemperature == "transrotational")
    {
        controlT_ = transrotational;
    }
    else if (controllingTemperature == "vibrational")
    {
        controlT_ = vibrational;
    }
    else
    {
        Info  << "Foam::Reaction2<Reaction2Thermo>::Reaction2(const speciesTable& species, "
              << "const HashPtrTable<Reaction2Thermo>& thermoDatabase, const dictionary& dict)"
              << "Controlling temperature " << controllingTemperature  << "is invalid." << nl << nl
              << "Valid Controlling temperature types are: " 
              << " chargeExchange, dissociation, exchange, impactDissociation," << nl 
              << " impactIonisation, associativeIonisation, transrotational, vibrational."
              << exit(FatalIOError);
    }
    // END BRAND NEW VINCENT **************************************************
    setLRhs
    (
        IStringStream(dict.lookup("reaction"))(),
        species_,
        lhs_,
        rhs_
    );
    
    setThermo(thermoDatabase);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Reaction2Thermo>
Foam::autoPtr<Foam::Reaction2<Reaction2Thermo> >
Foam::Reaction2<Reaction2Thermo>::New
(
    const speciesTable& species,
    const HashPtrTable<Reaction2Thermo>& thermoDatabase,
    Istream& is
)
{
    if (is.eof())
    {
        FatalIOErrorIn
        (
            "Reaction2<Reaction2Thermo>::New(const speciesTable&, "
            " const HashPtrTable<Reaction2Thermo>&, Istream&)",
            is
        )   << "Reaction2 type not specified" << nl << nl
            << "Valid Reaction2 types are :" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word reactionTypeName(is);

    typename IstreamConstructorTable::iterator cstrIter
        = IstreamConstructorTablePtr_->find(reactionTypeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "Reaction2<Reaction2Thermo>::New(const speciesTable&, "
            " const HashPtrTable<Reaction2Thermo>&, Istream&)",
            is
        )   << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<Reaction2<Reaction2Thermo> >
    (
        cstrIter()(species, thermoDatabase, is)
    );
}


template<class Reaction2Thermo>
Foam::autoPtr<Foam::Reaction2<Reaction2Thermo> >
Foam::Reaction2<Reaction2Thermo>::New
(
    const speciesTable& species,
    const HashPtrTable<Reaction2Thermo>& thermoDatabase,
    const dictionary& dict
)
{
    const word& reactionTypeName = dict.lookup("type");

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(reactionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Reaction2<Reaction2Thermo>::New"
            "("
                "const speciesTable&, "
                "const HashPtrTable<Reaction2Thermo>&, "
                "const dictionary&"
            ")"
        )   << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Reaction2<Reaction2Thermo> >
    (
        cstrIter()(species, thermoDatabase, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Reaction2Thermo>
void Foam::Reaction2<Reaction2Thermo>::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeKeyword("reaction") << reactionStr(reaction)
        << token::END_STATEMENT << nl;
}


template<class Reaction2Thermo>
Foam::scalar Foam::Reaction2<Reaction2Thermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class Reaction2Thermo>
Foam::scalar Foam::Reaction2<Reaction2Thermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class Reaction2Thermo>
Foam::scalar Foam::Reaction2<Reaction2Thermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class Reaction2Thermo>
const Foam::speciesTable& Foam::Reaction2<Reaction2Thermo>::species() const
{
    return species_;
}


template<class Reaction2Thermo>
const Foam::speciesTable& Foam::Reaction2<Reaction2Thermo>::gasSpecies() const
{
    notImplemented
    (
        "const speciesTable& gasSpecies() const"
        " for this reaction"
    );
    return *reinterpret_cast<speciesTable*>(0);
}


template<class Reaction2Thermo>
const Foam::List<typename Foam::Reaction2<Reaction2Thermo>::specieCoeffs>&
Foam::Reaction2<Reaction2Thermo>::glhs() const
{
    notImplemented
    (
        "inline const List<typename Reaction2<Reaction2Thermo>::specieCoeffs>&"
        "Reaction2<Reaction2Thermo>::glhs()"
    );
    return *reinterpret_cast<List<specieCoeffs>*>(0);
}


template<class Reaction2Thermo>
const Foam::List<typename Foam::Reaction2<Reaction2Thermo>::specieCoeffs>&
Foam::Reaction2<Reaction2Thermo>::grhs() const
{
    notImplemented
    (
        "inline const List<typename Reaction2<Reaction2Thermo>::specieCoeffs>&"
        "Reaction2<Reaction2Thermo>::grhs()"
    );
    return *reinterpret_cast<List<specieCoeffs>*>(0);
}

// ************************************************************************* //
