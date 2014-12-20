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

#include "IOstream.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystem, 0);

    defineRunTimeSelectionTable(coordinateSystem, origAxisDir);
    defineRunTimeSelectionTable(coordinateSystem, origRotation);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from origin and two axes
Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& dir
)
:
    name_(name),
    origin_(origin),
    axis_(axis),
    dir_(dir)
{}


// Construct from origin and a coordinate rotation
Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    name_(name),
    origin_(origin),
    axis_(cr.axis()),
    dir_(cr.direction())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    dir_(dict.lookup("direction"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateSystem::~coordinateSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::coordinateSystem::write(Ostream& os) const
{
    os  << nl << type()
        << " origin: " << origin()
        << " axis: " << axis() << " direction: " << direction() << nl;
}


void Foam::coordinateSystem::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl
        << "    origin " << origin() << token::END_STATEMENT << nl
        << "    axis " << axis() << token::END_STATEMENT << nl
        << "    direction " << direction() << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coordinateSystem& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const coordinateSystem& p");
    return os;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
