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

Description

\*---------------------------------------------------------------------------*/

#include "particle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Istream
template<class particleType>
particle<particleType>::particle
(
    const Cloud<particleType>& cloud, 
    Istream& is
)
:
    cloud_(cloud)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> position_ >> celli_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&position_),
            sizeof(position_) + sizeof(celli_)
        );
    }

    // Check state of Istream
    is.check("particle<particleType>::particle(Istream&)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class particleType>
Ostream& operator<<(Ostream& os, const particle<particleType>& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os << p.position_
           << token::SPACE << p.celli_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.position_),
            sizeof(p.position_) + sizeof(p.celli_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const particle<particleType>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
