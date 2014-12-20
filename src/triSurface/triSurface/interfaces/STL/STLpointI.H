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

#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline STLpoint::STLpoint()
{}


//- Construct from components
inline STLpoint::STLpoint(float x, float y, float z)
:
    x_(x), y_(y), z_(z)
{}


//- Construct from point
inline STLpoint::STLpoint(const point& pt)
:
    x_(float(pt.x())), y_(float(pt.y())), z_(float(pt.z()))
{}


//- Construct from istream
inline STLpoint::STLpoint(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline float STLpoint::x() const
{
    return x_;
}

inline float STLpoint::y() const
{
    return y_;
}

inline float STLpoint::z() const
{
    return z_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline STLpoint::operator point() const
{
    return point(x_, y_, z_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Istream& operator>>(Istream& is, STLpoint& stlp)
{
    stlp.x_ = readScalar(is);
    stlp.y_ = readScalar(is);
    stlp.z_ = readScalar(is);

    return is;
}


inline Ostream& operator<<(Ostream& os, const STLpoint& stlp)
{
    os  << stlp.x_ << token::SPACE
        << stlp.y_ << token::SPACE
        << stlp.z_;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
