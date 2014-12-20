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

#include "error.H"

#include "edgeMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from file
edgeMesh::edgeMesh(const fileName& fname)
{
    IFstream is(fname);

    if (is.good())
    {
        is >> points_ >> edges_;
    }
    else
    {
        FatalErrorIn("edgeMesh::edgeMesh(const fileName& fname)")
            << "cannot open file " << fname
            << abort(FatalError);
    }
}


// construct from Istream
edgeMesh::edgeMesh(Istream& is)
:
    points_(is),
    edges_(is)
{
    // Check state of Istream
    is.check("edgeMesh::edgeMesh(Istream&)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const edgeMesh& em)
{
    os  << em.points_ << nl << em.edges_ << endl;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const edgeMesh&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
