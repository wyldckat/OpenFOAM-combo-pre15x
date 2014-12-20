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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline string::string()
{}


inline string::string(const std::string& stdStr)
:
    std::string(stdStr)
{}


// Copy character array
inline string::string(const char* str)
:
    std::string(str)
{}


// Construct from a given number of characters in a character array
inline string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


// Construct from a single character
inline string::string(const char c)
:
    std::string(1, c)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline string string::operator()(const size_type i, const size_type n) const
{
    return substr(i, n);
}


inline string string::operator()(const size_type n) const
{
    return substr(0, n);
}


inline string::hash::hash()
{}


inline string::size_type string::hash::operator()(const string& key) const
{
    register size_type hashVal = 0;

    for (string::const_iterator iter=key.begin(); iter!=key.end(); ++iter)
    {
        hashVal = hashVal<<1 ^ *iter;
    }

    return hashVal;
}


inline string::size_type string::hash::operator()
(
    const string& key,
    const size_type tableSize
) const
{
    return ::abs(operator()(key)) % tableSize;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
