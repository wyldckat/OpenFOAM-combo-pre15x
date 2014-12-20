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
    Reads a char from an input stream, for a given version
    number and File format. If an ascii File is being read, then the line
    numbers are counted and an erroneous read is reported.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "char.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Istream& operator>>(Istream& is, char& c)
{
    is.read(c);
    is.check("Istream& operator>>(Istream& is, char& c)");
    return is;
}


Ostream& operator<<(Ostream& os, const char c)
{
    os.write(c);
    os.check("Ostream& operator<<(Ostream&, const char)");
    return os;
}


char readChar(Istream& is)
{
   char c;
   is.read(c);
   return c;
}


Ostream& operator<<(Ostream& os, const char* s)
{
    os.write(s);
    os.check("Ostream& operator<<(Ostream&, const char*)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
