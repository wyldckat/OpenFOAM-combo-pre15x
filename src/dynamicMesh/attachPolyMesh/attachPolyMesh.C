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

\*---------------------------------------------------------------------------*/

#include "attachPolyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::attachPolyMesh::attachPolyMesh(const IOobject& io)
:
    polyMesh(io)
{}


// Construct from components without boundary.
// Boundary is added using addPatches() member function
Foam::attachPolyMesh::attachPolyMesh
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


//- Construct from cell shapes
Foam::attachPolyMesh::attachPolyMesh
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::attachPolyMesh::~attachPolyMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Attach mesh
void Foam::attachPolyMesh::attach(const bool removeEmptyPatches)
{
    if (polyMesh::debug || polyMesh::morphDebug)
    {
        Info<< "void attachPolyMesh::attach(): "
            << "Attaching mesh" << endl;
    }

    // Save current file instance
    const fileName oldInst = cellsInstance();

    // Reset morph time index
    setMorphTimeIndex(-1);

    // Execute all polyMeshModifiers
    updateTopology();

    const pointField p = oldAllPoints();

    movePoints(p);

    if (polyMesh::debug || polyMesh::morphDebug)
    {
        Info << "Clearing mesh." << endl;
    }

    removePointZones();
    removeFaceZones();
    removeCellZones();
    removeMeshModifiers();

    if (removeEmptyPatches)
    {
        // Re-do the boundary patches, removing the ones with zero size
        const polyBoundaryMesh& oldPatches = boundaryMesh();

        List<polyPatch*> newPatches(oldPatches.size());
        label nNewPatches = 0;

        forAll (oldPatches, patchI)
        {
            if (oldPatches[patchI].size() > 0)
            {
                newPatches[nNewPatches] =
                    oldPatches[patchI].clone
                    (
                        boundaryMesh(),
                        nNewPatches,
                        oldPatches[patchI].size(),
                        oldPatches[patchI].start()
                    ).ptr();

                nNewPatches++;
            }
            else
            {
                if (polyMesh::debug || polyMesh::morphDebug)
                {
                    Info<< "Removing zero-sized patch " << patchI
                        << " named " << oldPatches[patchI].name() << endl;
                }
            }
        }

        newPatches.setSize(nNewPatches);

        removeBoundary();
        addPatches(newPatches);
    }

    // Reset the file instance to overwrite the original mesh
    setInstance(oldInst);

    if (polyMesh::debug || polyMesh::morphDebug)
    {
        Info<< "void attachPolyMesh::attach(): "
            << "Finished attaching mesh" << endl;
    }

    checkMesh();
}


// Write after removing files from old mesh
bool Foam::attachPolyMesh::write
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (polyMesh::debug)
    {
        Info<< "void attachPolyMesh::attach(): "
            << "Removing old mesh files from "
            << cellsInstance() << endl;
    }

    removeFiles(cellsInstance());

    if (polyMesh::debug)
    {
        Info<< "Writing attachPolyMesh to " << cellsInstance() << endl;
    }

    return polyMesh::write(fmt, ver, cmp);
}


bool Foam::attachPolyMesh::write() const
{
    return regIOobject::write();
}


// ************************************************************************* //
