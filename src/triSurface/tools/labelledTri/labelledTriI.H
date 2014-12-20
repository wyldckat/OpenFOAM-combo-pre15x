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

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline labelledTri::labelledTri()
:
    region_(-1)
{}


//- Construct from components
inline labelledTri::labelledTri
(
    const label A,
    const label B,
    const label C,
    const label region
)
:
    triFace(A, B, C),
    region_(region)
{}


inline labelledTri::labelledTri(Istream& is)
{
    // Read beginning of labelledTri point pair
    is.readBegin("labelledTri");

    is >> (triFace&)(*this) >> region_;

    // Read end of labelledTri point pair
    is.readEnd("labelledTri");

    // Check state of Istream
    is.check("labelledTri::labelledTri(Istream& is)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label labelledTri::region() const
{
    return region_;
}

inline label& labelledTri::region()
{
    return region_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Istream& operator>>(Istream& is, labelledTri& t)
{
    // Read beginning of labelledTri point pair
    is.readBegin("labelledTri");

    is  >> (triFace&)t >> t.region_;

    // Read end of labelledTri point pair
    is.readEnd("labelledTri");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, labelledTri&)");

    return is;
}


inline Ostream& operator<<(Ostream& os, const labelledTri& t)
{
    os  << token::BEGIN_LIST
        << (const triFace&)t << token::SPACE << t.region_
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
