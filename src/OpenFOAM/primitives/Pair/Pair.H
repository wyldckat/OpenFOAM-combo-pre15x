/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    Pair

Description
    An ordered Pair.


\*---------------------------------------------------------------------------*/

#ifndef Pair_H
#define Pair_H

#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class Pair;

template<class Type>
Istream& operator>>(Istream&, Pair<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const Pair<Type>&);

/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Pair
{
    // Private data

       Type first_, second_;


public:

    // Constructors

        //- Null constructor for lists
        inline Pair()
        {}

        //- Construct from components
        inline Pair(const Type& first, const Type& second)
        :
            first_(first),
            second_(second)
        {}

        //- Construct from Istream
        inline Pair(Istream& is)
        {
            // Read beginning of pair
            is.readBegin("pair");

            is >> first_ >> second_;

            // Read end of pair
            is.readEnd("pair");

            // Check state of Istream
            is.check("Pair::Pair(Istream&)");
        }


    // Member Functions

        //- Return first
        inline Type first() const
        {
            return first_;
        }

        //- Return first
        inline Type& first()
        {
            return first_;
        }

        //- Return second
        inline Type second() const
        {
            return second_;
        }

        //- Return second
        inline Type& second()
        {
            return second_;
        }

        //- Return reverse pair
        inline Pair<Type> reversePair() const
        {
            return Pair<Type>(second_, first_);
        }


    // Friend Operators

        inline friend bool operator==(const Pair<Type>& a, const Pair<Type>& b)
        {
            return
            (
                (a.first_ == b.first_) && (a.second_ == b.second_)
            );
        }

        inline friend bool operator!=(const Pair<Type>& a, const Pair<Type>& b)
        {
            return (!(a == b));
        }


    // IOstream Operators

        friend Istream& operator>> <Type>(Istream& is, Pair<Type>& p);
        friend Ostream& operator<< <Type>(Ostream& os, const Pair<Type>& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Istream& operator>>(Istream& is, Pair<Type>& p)
{
    // Read beginning of Pair<Type>
    is.readBegin("Pair<Type>");

    is >> p.first_ >> p.second_;

    // Read end of Pair<Type>
    is.readEnd("Pair<Type>");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, Pair<Type>&)");

    return is;
}

template<class Type>
Ostream& operator<<(Ostream& os, const Pair<Type>& p)
{
    os  << token::BEGIN_LIST
        << p.first_ << token::SPACE
        << p.second_
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const Pair<Type>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
