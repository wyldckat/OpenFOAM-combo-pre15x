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
    A class for handling file names.
    A fileName can be  constructed from a char* or a word.
    FileNames may be concaternated by adding a '/' separator.
    FileNames may be decomposed into the path, name or component
    list. FileNames may be interogated for type and access mode.

\*---------------------------------------------------------------------------*/

#include "fileName.H"
#include "wordList.H"
#include "debug.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int fileName::debug(debug::debugSwitch("fileName", 0));
const fileName fileName::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fileName::fileName(const wordList& wrdList)
{
    if (wrdList.size() != 0)
    {
        forAll(wrdList, i)
        {
            operator=((*this)/wrdList[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return file name extension
word fileName::name() const
{
    size_type i = rfind('/');

    if (i > 0)
    {
        return substr(i+1, size()-i);
    }
    else
    {
        return word::null;
    }
}


// Return file name extension
fileName fileName::path() const
{
    size_type i = rfind('/');

    if (i > 0)
    {
        return substr(0, i);
    }
    else
    {
        return fileName::null;
    }
}


// Return file name less extension (part before last .)
fileName fileName::lessExt() const
{
    size_type i = find_last_of("./");

    if (i <= 0 || operator[](i) == '/')
    {
        return *this;
    }
    else
    {
        return substr(0, i);
    }
}


// Return file name extension
word fileName::ext() const
{
    size_type i = find_last_of("./");

    if (i <= 0 || operator[](i) == '/')
    {
        return word::null;
    }
    else
    {
        return substr(i+1, size()-i);
    }
}


// Return the components of the file name as a wordList
wordList fileName::components(const char delimiter) const
{
    wordList wrdList(20);

    size_type start=0, end=0;
    label nWords=0;

    while ((end = find(delimiter, start)) != npos)
    {
        wrdList[nWords++] = substr(start, end-start);

        if (nWords == wrdList.size())
        {
            wrdList.setSize(2*wrdList.size());
        }

        start = end+1;
    }

    wrdList[nWords++] = substr(start, size()-start);

    wrdList.setSize(nWords);

    return wrdList;
}


// Return a component of the file name
word fileName::component(const size_type cmpt, const char delimiter) const
{
    return components(delimiter)[cmpt];
}


fileName::Type fileName::type() const
{
    return ::Foam::type(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void fileName::operator=(const fileName& q)
{
    string::operator=(q);
}


void fileName::operator=(const word& q)
{
    string::operator=(q);
}


void fileName::operator=(const string& q)
{
    string::operator=(q);
    stripInvalid();
}


void fileName::operator=(const std::string& q)
{
    string::operator=(q);
    stripInvalid();
}


void fileName::operator=(const char* q)
{
    string::operator=(q);
    stripInvalid();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

fileName operator/(const string& a, const string& b)
{
    if (a.size() > 0)       // First string non-null
    {
        if (b.size() > 0)   // Second string non-null
        {
            return fileName(a + '/' + b);
        }
        else                // Second string null
        {
            return a;
        }
    }
    else                    // First string null
    {
        if (b.size() > 0)   // Second string non-null
        {
            return b;
        }
        else                // Second string null
        {
            return fileName();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
