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

\*---------------------------------------------------------------------------*/

#include "morphMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::morphMesh::morphMesh(const IOobject& io)
:
    polyMesh(io)
{}


// Construct from components without boundary.
// Boundary is added using addPatches() member function
Foam::morphMesh::morphMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    polyMesh
    (
        io,
        points,
        faces,
        cells
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::morphMesh::~morphMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Add boundary patches
void Foam::morphMesh::removeBoundary()
{
    polyMesh::removeBoundary();
}


// Add boundary patches
void Foam::morphMesh::addPatches(const List<polyPatch*>& patches)
{
    polyMesh::addPatches(patches);
}


// ************************************************************************* //
