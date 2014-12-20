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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Pair

Description
    An ordered Pair.

\*---------------------------------------------------------------------------*/

#ifndef Pair_H
#define Pair_H

#include "FixedList.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Pair
:
    public FixedList<Type, 2>
{

public:

    // Constructors

        //- Null constructor for lists
        inline Pair()
        {}

        //- Construct from components
        inline Pair(const Type& f, const Type& s)
        {
            first() = f;
            second() = s;
        }

        //- Construct from Istream
        inline Pair(Istream& is)
        :
            FixedList<Type, 2>(is)
        {}


    // Member Functions

        //- Return first
        inline Type first() const
        {
            return this->operator[](0);
        }

        //- Return first
        inline Type& first()
        {
            return this->operator[](0);
        }

        //- Return second
        inline Type second() const
        {
            return this->operator[](1);
        }

        //- Return second
        inline Type& second()
        {
            return this->operator[](1);
        }

        //- Return reverse pair
        inline Pair<Type> reversePair() const
        {
            return Pair<Type>(second(), first());
        }


    // Friend Operators

        inline friend bool operator==(const Pair<Type>& a, const Pair<Type>& b)
        {
            return
            (
                (a.first() == b.first()) && (a.second() == b.second())
            );
        }

        inline friend bool operator!=(const Pair<Type>& a, const Pair<Type>& b)
        {
            return (!(a == b));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
