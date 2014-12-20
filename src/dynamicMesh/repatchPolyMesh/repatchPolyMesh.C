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
    A mesh which allows changes in the patch distribution of the
    faces.  The change in patching is set using changePatchID.  For a
    boundary face, a new patch ID is given.  If the face is internal,
    it will be added to the first patch and its opposite to the second
    patch (take care with face orientation!).

\*---------------------------------------------------------------------------*/

#include "repatchPolyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(repatchPolyMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::repatchPolyMesh::repatchPolyMesh(const IOobject& io)
:
    polyMesh(io),
    meshMod_(*this)
{}


// Construct from components
Foam::repatchPolyMesh::repatchPolyMesh
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
    ),
    meshMod_(*this)
{}


// Construct from cell shapes
Foam::repatchPolyMesh::repatchPolyMesh
(
    const IOobject& io,
    const pointField& points,
    const cellShapeList& shapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const wordList& boundaryPatchTypes,
    const word& defaultBoundaryPatchType,
    const wordList& boundaryPatchPhysicalTypes
)
:
    polyMesh
    (
        io,
        points,
        shapes,
        boundaryFaces,
        boundaryPatchNames,
        boundaryPatchTypes,
        defaultBoundaryPatchType,
        boundaryPatchPhysicalTypes
    ),
    meshMod_(*this)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::repatchPolyMesh::changePatchID
(
    const label faceID,
    const label patchID
)
{
    if (debug)
    {
        // Check that the request is possible
        if
        (
            faceID >= allFaces().size()
         || patchID >= boundaryMesh().size()
         || !isInternalFace(faceID)
        )
        {
            FatalErrorIn
            (
                "void Foam::repatchPolyMesh::changePatchID\n"
                "(\n"
                "    const label faceID,\n"
                "    const label patchID\n"
                ")\n"
            )   << "Error in definition.  faceID: " << faceID
                << " patchID: " << patchID << ".  "
                << "Labels out of range or internal face."
                << abort(FatalError);
        }
    }

    const label zoneID = faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod_.setAction
    (
        polyModifyFace
        (
            allFaces()[faceID],                 // face
            faceID,                             // face ID
            allOwner()[faceID],                 // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


void Foam::repatchPolyMesh::setFaceZone
(
    const label faceID,
    const label zoneID,
    const bool zoneFlip
)
{
    if (debug)
    {
        // Check that the request is possible
        if (faceID > allFaces().size())
        {
            FatalErrorIn
            (
                "void Foam::repatchPolyMesh::setFaceZone"
                "(\n"
                "    const label faceID,\n"
                "    const label zoneID,\n"
                "    const bool flip\n"
                ")\n"
            )   << "Error in definition.  faceID: " << faceID
                << "out of range."
                << abort(FatalError);
        }
    }

    meshMod_.setAction
    (
        polyModifyFace
        (
            allFaces()[faceID],                 // face
            faceID,                             // face ID
            allOwner()[faceID],                 // owner
            allNeighbour()[faceID],             // neighbour
            false,                              // flip flux
            boundaryMesh().whichPatch(faceID),  // patch ID
            true,                               // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


void Foam::repatchPolyMesh::changeAnchorPoint
(
    const label faceID,
    const label fp
)
{
    if (debug)
    {
        // Check that the request is possible
        if (faceID > allFaces().size())
        {
            FatalErrorIn
            (
                "void Foam::repatchPolyMesh::setFaceZone"
                "(\n"
                "    const label faceID,\n"
                "    const label zoneID,\n"
                "    const bool flip\n"
                ")\n"
            )   << "Error in definition.  faceID: " << faceID
                << "out of range."
                << abort(FatalError);
        }
    }

    const face& f = allFaces()[faceID];

    if ((fp < 0) || (fp >= f.size()))
    {
        FatalErrorIn
        (
            "void Foam::repatchPolyMesh::changeAnchorPoint"
            "(\n"
            "    const label faceID,\n"
            "    const label fp\n"
            ")\n"
        )   << "Error in definition.  Face point: " << fp
            << "indexes out of face " << f
            << abort(FatalError);
    }

    label patchID = boundaryMesh().whichPatch(faceID);

    const label zoneID = faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    if (fp == 0)
    {
        // Do dummy modify to keep patch ordering.
        meshMod_.setAction
        (
            polyModifyFace
            (
                f,                                  // face
                faceID,                             // face ID
                allOwner()[faceID],                 // owner
                -1,                                 // neighbour
                false,                              // flip flux
                patchID,                            // patch ID
                false,                              // remove from zone
                zoneID,                             // zone ID
                zoneFlip                            // zone flip
            )
        );
    }
    else
    {
        // Construct new face with fp as first point.

        face newFace(f.size());

        label fVert = fp;

        for (label i = 0; i < f.size(); i++)
        {
            newFace[i] = f[fVert++];

            if (fVert == f.size())
            {
                fVert = 0;
            }
        }


        meshMod_.setAction
        (
            polyModifyFace
            (
                newFace,                            // face
                faceID,                             // face ID
                allOwner()[faceID],                 // owner
                -1,                                 // neighbour
                false,                              // flip flux
                patchID,                            // patch ID
                false,                              // remove from zone
                zoneID,                             // zone ID
                zoneFlip                            // zone flip
            )
        );
    }
}


void Foam::repatchPolyMesh::repatch()
{
    setMorphTimeIndex(time().timeIndex());
    updateTopology(meshMod_);

    // Clear topo change for the next operation
    meshMod_.clear();
}


// ************************************************************************* //
