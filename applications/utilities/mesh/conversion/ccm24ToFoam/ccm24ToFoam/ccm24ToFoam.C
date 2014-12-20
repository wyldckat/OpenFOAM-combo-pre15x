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
    Reads CCM/NGeom files as written by Prostar/ccm.

    - does polyhedral mesh
    - loses patch names (creates 'patch0', 'patch1' etc. instead)
    - does not handle 'interfaces' (couples)
    - does not handle cyclics
    - does not do data

    Uses libccmiosrc/test/Mesh* classes to do actual reading. At our level we
    only add the compacting of vertices and cells (vertices might be
    non-compact, cells might start from 1 but probably are compact?)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "SortableList.H"

#include "Mesh.h"
#include "MeshCCMIO.h"
#include "MeshNGEOM.h"
#include "PostCCMIO.h"
#include "PostNPOST.h"
#include "IOException.h"
#include <iostream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Determine order for faces:
// - upper-triangular order for internal faces
// - external faces after internal faces (not really nessecary since already
//   so)
labelList getFaceOrder
(
    const cellList& cells,
    const labelList& owner,
    const labelList& neighbour,
    const label nInternalFaces
)
{
    labelList oldToNew(owner.size(), -1);

    // First unassigned internal face
    label internalFaceI = 0;

    // First unassigned boundary face
    label boundaryFaceI = nInternalFaces;

    forAll(cells, cellI)
    {
        const labelList& cFaces = cells[cellI];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            label nbrCellI = neighbour[faceI];

            if (nbrCellI != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCellI == cellI)
                {
                    nbrCellI = owner[faceI];
                }

                if (cellI < nbrCellI)
                {
                    // CellI is master
                    nbr[i] = nbrCellI;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Assign to first free boundary face spot.
                nbr[i] = -1;

                oldToNew[faceI] = boundaryFaceI++;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNew[cFaces[nbr.indices()[i]]] = internalFaceI++;
            }
        }
    }


    // Check done all faces.
    forAll(oldToNew, faceI)
    {
        if (oldToNew[faceI] == -1)
        {
            FatalErrorIn("getFaceOrder") << "Did not determine new position"
                << " for face " << faceI
                << abort(FatalError);
        }
    }


    return oldToNew;
}


// Use the map from old to new numbering to update the face information.
void reorderFaces
(
    const labelList& oldToNew,
    cellList& cells,
    faceList& faces,
    labelList& owner,
    labelList& neighbour
)
{
    // Renumber faces in cells
    forAll(cells, cellI)
    {
        cell& cFaces = cells[cellI];

        forAll(cFaces, i)
        {
            cFaces[i] = oldToNew[cFaces[i]];
        }
    }

    // Reorder faces.
    faceList oldFaces(faces);
    labelList oldOwner(owner);
    labelList oldNeighbour(neighbour);

    forAll(oldToNew, faceI)
    {
        label newFaceI = oldToNew[faceI];

        faces[newFaceI].transfer(oldFaces[faceI]);
        owner[newFaceI] = owner[faceI];
        neighbour[newFaceI] = neighbour[faceI];
    }
}



// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("ngeom/ccm file");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName ccmFile(args.args()[3]);

    if (!exists(ccmFile))
    {
        FatalErrorIn(args.executable())
            << "Cannot read file " << ccmFile
            << exit(FatalError);
    }

    word ccmExt = ccmFile.ext();

    autoPtr<Mesh> meshPtr(NULL);

    if (ccmExt == "ccm" || ccmExt == "sff")
    {
        meshPtr.reset(new MeshCCMIO(ccmFile.c_str()));
    }
    else if (ccmExt == "ngeom")
    {
        meshPtr.reset(new MeshNGEOM(ccmFile.c_str()));
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Illegal extension " << ccmExt << " for file " << ccmFile << nl
            << "Allowed extensions are '.ccm', '.sff' or '.ngeom'"
            << exit(FatalError);
    }

    const Mesh& ccMesh = meshPtr();

    //
    // Convert MeshVertex to pointField
    //

    const Mesh::VertexArray& ccVerts = ccMesh.vertices;

    // Vertex 0 never used so store with one offset.

    Info<< "Number of (uncompacted) points:" << label(ccVerts.size()) << endl;

    pointField points(ccVerts.size());

    forAll(ccVerts, pointI)
    {
        const ::MeshVertex& ccVert = ccVerts[pointI];

        points[pointI] =
            Foam::point
            (
                ccVert.coord[0],
                ccVert.coord[1],
                ccVert.coord[2]
            );
    }


    //
    // MeshFace to faceList
    //

    const Mesh::FaceArray& ccIntFaces = ccMesh.internalFaces;
    const Mesh::FaceArray& ccExtFaces = ccMesh.boundaryFaces;
    const Mesh::CellArray& ccCells = ccMesh.cells;

    Info<< "Number of internal faces:" << label(ccIntFaces.size()) << endl;
    Info<< "Number of boundary faces:" << label(ccExtFaces.size()) << endl;
    Info<< "Number of (uncompacted) cells:" << label(ccCells.size()) << endl;

    // Foam storage
    faceList faces(ccIntFaces.size() + ccExtFaces.size());
    labelList owner(faces.size(), -1);
    labelList neighbour(faces.size(), -1);
    labelList region(faces.size(), -1);

    label faceI = 0;

    // Internal faces
    for
    (
        Mesh::FaceArray::const_iterator iter = ccIntFaces.begin();
        iter != ccIntFaces.end();
        ++iter
    )
    {
        const ::MeshFace& ccFace = *iter;

        face& f = faces[faceI];

        f.setSize(ccFace.nVerts);

        forAll(f, fp)
        {
            f[fp] = ccFace.verts[fp];
        }

        owner[faceI] = ccFace.cells[0];
        neighbour[faceI] = ccFace.cells[1];

        if (neighbour[faceI] < 0)
        {
            FatalErrorIn(args.executable())
                << "Internal face has invalid neighbour cell." << nl
                << "Face:" << faceI << " verts:" << f
                << " owner:" << owner[faceI]
                << " neighbour:" << neighbour[faceI]
                << abort(FatalError);
        }

        faceI++;
    }

    // External faces
    label maxPatch = labelMin;

    for
    (
        Mesh::FaceArray::const_iterator iter = ccExtFaces.begin();
        iter != ccExtFaces.end();
        ++iter
    )
    {
        const ::MeshFace& ccFace = *iter;

        face& f = faces[faceI];

        f.setSize(ccFace.nVerts);

        forAll(f, fp)
        {
            f[fp] = ccFace.verts[fp];
        }

        owner[faceI] = ccFace.cells[0];

        region[faceI] = ccFace.boundary;

        maxPatch = max(maxPatch, region[faceI]);

        if (ccFace.cells[1] >= 0)
        {
            FatalErrorIn(args.executable())
                << "Boundary face has valid neighbour cell." << nl
                << "Face:" << faceI << " verts:" << f
                << " region:" << region[faceI]
                << " owner:" << owner[faceI]
                << " neighbour:" << ccFace.cells[1]
                << abort(FatalError);
        }

        faceI++;
    }

    const label nPatches = maxPatch+1;

    Info<< "Number of patches:" << nPatches << endl;


    // We now have extracted all info from CCMIO:
    // - coordinates (points)
    // - face to point addressing (faces)
    // - face to cell addressing (owner, neighbour, region)
    //
    // Problem is that points and cells are non-compact so compact


    // Compact points.
    {
        labelList pointMap(points.size(), -1);

        pointField newPoints(points.size());
        label newPointI = 0;

        forAll(faces, faceI)
        {
            face& f = faces[faceI];

            forAll(f, fp)
            {
                label pointI = f[fp];

                if (pointMap[pointI] == -1)
                {
                    pointMap[pointI] = newPointI;

                    newPoints[newPointI] = points[pointI];

                    f[fp] = newPointI;

                    newPointI++;
                }
                else
                {
                    f[fp] = pointMap[pointI];
                }
            }
        }

        Info<< "Compacted points from " << points.size() << " down to "
            << newPointI << endl;

        newPoints.setSize(newPointI);

        points.transfer(newPoints);
    }


    // Compact cells
    label nCells = 0;
    {
        labelList cellMap(owner.size(), -1);

        forAll(owner, faceI)
        {
            label own = owner[faceI];

            if (cellMap[own] == -1)
            {
                cellMap[own] = nCells++;
            }
            owner[faceI] = cellMap[own];


            label nei = neighbour[faceI];

            if (nei >= 0)
            {
                label nei = neighbour[faceI];

                if (cellMap[nei] == -1)
                {
                    cellMap[nei] = nCells++;
                }
                neighbour[faceI] = cellMap[nei];
            }
        }
        Info<< "Compacted cells from " << label(ccCells.size()) << " down to "
            << nCells << endl;
    }



    // Create cells (inverse of face-to-cell addressing)
    cellList cells(nCells);

    {
        // Create cell-faces addressing.
        List<DynamicList<label> > cellFaces(cells.size());

        forAll(owner, faceI)
        {
            cellFaces[owner[faceI]].append(faceI);

            if (neighbour[faceI] >= 0)
            {
                cellFaces[neighbour[faceI]].append(faceI);
            }
        }

        // Transfer cellFaces to proper cells
        forAll(cellFaces, cellI)
        {
            cellFaces[cellI].shrink();
            cells[cellI].transfer(cellFaces[cellI]);
            cellFaces[cellI].clear();
        }
    }

    //
    // Basic mesh info complete. Now convert to Foam convention:
    // - owner < neighbour
    // - face vertices such that normal points away from owner
    // - order faces: upper-triangular for internal faces; boundary faces after
    //   internal faces
    //


    // Set owner/neighbour so owner < neighbour
    forAll(neighbour, faceI)
    {
        label nbr = neighbour[faceI];
        label own = owner[faceI];

        if (nbr >= cells.size() || own < 0 || own >= cells.size())
        {
            FatalErrorIn(args.executable())
                << "face:" << faceI 
                << " nbr:" << nbr
                << " own:" << own
                << abort(FatalError);
        }

        if (nbr >= 0)
        {
            if (nbr < own)
            {
                owner[faceI] = neighbour[faceI];
                neighbour[faceI] = own;
                reverse(faces[faceI]);
            }
        }
    }

    // Determine face order for upper-triangular ordering and internal/external
    // face ordering
    labelList oldToNew(getFaceOrder(cells, owner, neighbour, ccIntFaces.size()));

    // Get starting face label of patches
    labelList patchSizes(nPatches, 0);

    forAll(region, faceI)
    {
        if (region[faceI] >= 0)
        {
            patchSizes[region[faceI]]++;
        }
    }
    Info<< "patchSizes:" << patchSizes << endl;


    labelList start(nPatches, 0);

    label meshFaceI = ccIntFaces.size();

    for (label patchI = 0; patchI < nPatches; patchI++)
    {
        start[patchI] = meshFaceI;

        meshFaceI += patchSizes[patchI];
    }
    Info<< "patches starting at faces:" << start << endl;

    forAll(region, faceI)
    {
        if (region[faceI] >= 0)
        {
            oldToNew[faceI] = start[region[faceI]]++;
        }
    }


    // Reorder faces accordingly
    reorderFaces(oldToNew, cells, faces, owner, neighbour);

    // Construct polyMesh (without patches)
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        points,
        faces,
        cells
    );

    // Create single empty patch
    List<polyPatch*> newPatches(nPatches);

    {
        label meshFaceI = ccIntFaces.size();

        forAll(newPatches, patchI)
        {
            newPatches[patchI] =
                new emptyPolyPatch
                (
                    "patch" + Foam::name(patchI),
                    start[patchI] - meshFaceI,
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh()
                );

            meshFaceI = start[patchI];
        }
    }

    mesh.addPatches(newPatches);

    Info<< "Writing polyMesh to " << runTime.constant()/polyMesh::defaultRegion
        << "..." << nl << endl;

    mesh.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
