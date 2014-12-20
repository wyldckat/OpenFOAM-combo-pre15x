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

#include "morphFvMesh.H"
#include "parallelInfo.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Read field
template<class GeoField>
void Foam::morphFvMesh::readFields
(
    const IOobjectList& objects,
    ptrList<GeoField>& fields
) const
{
    // Search list of objects for volScalarFields
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Construct the vol scalar fields
    fields.setSize(fieldObjects.size());

    for
    (
        IOobjectList::iterator iter = fieldObjects.begin();
        iter != fieldObjects.end();
        ++iter
    )
    {
        fields.hook(new GeoField(*iter(), *this));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::morphFvMesh::morphFvMesh(const IOobject& io)
:
    fvMesh(io)
{}


// Construct from components without boundary.
// Boundary is added using addPatches() member function
Foam::morphFvMesh::morphFvMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    fvMesh(io, points, faces, cells)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::morphFvMesh::~morphFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Add boundary patches
void Foam::morphFvMesh::addFvPatches(const List<polyPatch*>& patches)
{
    fvMesh::addFvPatches(patches);
}


//- Remove fvPatches
void Foam::morphFvMesh::removeFvBoundary()
{
    fvMesh::removeFvBoundary();
}


void Foam::morphFvMesh::updateTopology(const polyTopoChange& changer)
{
    // Like fvMesh::updateTopology but with explicitly provided changes.
    polyMesh::updateTopology(changer);

    handleMorph();
}


bool  Foam::morphFvMesh::setMorphTimeIndex(const label newTimeIndex) const
{
    return polyMesh::setMorphTimeIndex(newTimeIndex);
}


// ************************************************************************* //
