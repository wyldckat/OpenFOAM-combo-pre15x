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
    Virtual base class for mesh modifiers.

\*---------------------------------------------------------------------------*/

#include "polyMeshModifier.H"
#include "dictionary.H"
#include "polyMeshMorphEngine.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshModifier, 0);

    defineRunTimeSelectionTable(polyMeshModifier, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyMeshModifier::polyMeshModifier
(
    const word& name,
    const label index,
    const polyMeshMorphEngine& mme,
    const bool act
)
:
    name_(name),
    index_(index),
    morphEngine_(mme),
    active_(act)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMeshModifier::~polyMeshModifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMeshMorphEngine& Foam::polyMeshModifier::morphEngine() const
{
    return morphEngine_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyMeshModifier& pmm)
{
    pmm.write(os);
    os.check("Ostream& operator<<(Ostream& f, const polyMeshModifier& pmm)");
    return os;
}


// ************************************************************************* //
