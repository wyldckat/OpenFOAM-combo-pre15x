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
    Reads CCM/NGeom files as written by Prostar/ccm.

    - does polyhedral mesh
    - does not handle 'interfaces' (couples)
    - does not handle cyclics
    - does not do data
    - not tested for patch names (these seem to go as 'Label' into adf file?)

    Uses libccmiosrc/test/Mesh* classes to do actual reading. At our level we
    only add the compacting of vertices and cells (vertices might be
    non-compact, cells might start from 1 but probably are compact?)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "SortableList.H"

#include "Mesh.h"
#include "MeshCCMIO.h"
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
    const labelList& region,
    const label nPatches,
    labelList& patchSizes
)
{
    labelList oldToNew(owner.size(), -1);

    // First unassigned face
    label newFaceI = 0;

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
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNew[cFaces[nbr.indices()[i]]] = newFaceI++;
            }
        }
    }


    // Pick up all patch faces in patch face order. (note: loops over all
    // faces for all patches. Not very efficient)
    patchSizes.setSize(nPatches);
    patchSizes = 0;

    for (label patchI = 0; patchI < nPatches; patchI++)
    {
        forAll(region, faceI)
        {
            if (region[faceI] == patchI)
            {
                oldToNew[faceI] = newFaceI++;
                patchSizes[patchI]++;
            }
        }
    }

    // Check done all faces.
    forAll(oldToNew, faceI)
    {
        if (oldToNew[faceI] == -1)
        {
            FatalErrorIn("getFaceOrder") << "Did not determine new position"
                << " for face " << faceI << endl
                << "Face owner cell: " << owner[faceI]
                << ", face neighbour cell: " << neighbour[faceI] << endl
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

    // Foam storage of mesh data
    faceList faces;
    labelList owner;
    labelList neighbour;
    labelList region;
    labelList cellId;
    pointField points;
    List<word> patchNames;
    List<word> patchTypes;
    label nInternalFaces = -1;
    
    {
        fileName ccmFile(args.args()[3]);

        if (!exists(ccmFile))
        {
            FatalErrorIn(args.executable())
                << "Cannot read file " << ccmFile
                << exit(FatalError);
        }

        word ccmExt = ccmFile.ext();

        if (ccmExt != "ccm" && ccmExt != "sff")
        {
            FatalErrorIn(args.executable())
                << "Illegal extension " << ccmExt << " for file " << ccmFile
                << nl << "Allowed extensions are '.ccm', '.sff'"
                << exit(FatalError);
        }

        MeshCCMIO ccMesh(ccmFile.c_str());


        //
        // Convert MeshVertex to pointField
        //

        const Mesh::VertexArray& ccVerts = ccMesh.vertices;

        Info<< "Number of (uncompacted) points:" << label(ccVerts.size())
            << endl;

        points.setSize(ccVerts.size());

        for (unsigned int pointI = 0; pointI < ccVerts.size(); pointI++)
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
        Info<< "Number of (uncompacted) cells:" << label(ccCells.size())
            << endl;

        // Size Foam storage
        nInternalFaces = ccIntFaces.size();
        faces.setSize(nInternalFaces + ccExtFaces.size());
        owner.setSize(faces.size());
        owner = -1;
        neighbour.setSize(faces.size());
        neighbour = -1;
        region.setSize(faces.size());
        region = -1;
        cellId.setSize(ccCells.size());
        cellId = -1;
        
        // Cells ids
        label cellI = 0;
        for
        (
            Mesh::CellArray::const_iterator iter = ccCells.begin();
            iter != ccCells.end();
            ++iter, cellI++
        )
        {
            const ::MeshCell& ccCell = *iter;

            cellId[cellI] = ccCell.id;
        }


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

            f.setSize(ccFace.verts.size());

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

        // Map from ccm boundary id to Foam patch
        Map<label> ccmToPatch(ccMesh.boundaryRegionInfo.size());

        // Patches
        patchNames.setSize(ccMesh.boundaryRegionInfo.size());
        patchTypes.setSize(patchNames.size());
        forAll(patchNames, pnI)
        {
            ccmToPatch.insert(ccMesh.boundaryRegionInfo[pnI].boundary, pnI);
            patchNames[pnI] = ccMesh.boundaryRegionInfo[pnI].label;
            patchTypes[pnI] = ccMesh.boundaryRegionInfo[pnI].boundaryType;

            Info<< "Patch:" << pnI
                << " ccm boundary id:"
                << ccMesh.boundaryRegionInfo[pnI].boundary
                << " name:" << patchNames[pnI]
                << " type:" << patchTypes[pnI] << endl;
        }


        // External faces. Remap region.

        for
        (
            Mesh::FaceArray::const_iterator iter = ccExtFaces.begin();
            iter != ccExtFaces.end();
            ++iter
        )
        {
            const ::MeshFace& ccFace = *iter;

            face& f = faces[faceI];

            f.setSize(ccFace.verts.size());

            forAll(f, fp)
            {
                f[fp] = ccFace.verts[fp];
            }

            owner[faceI] = ccFace.cells[0];

            region[faceI] = ccmToPatch[ccFace.boundary];

            if (ccFace.cells[1] >= 0)
            {
                FatalErrorIn(args.executable())
                    << "Boundary face has valid neighbour cell." << nl
                    << "Face:" << faceI << " verts:" << f
                    << " region:" << region[faceI]
                    << " owner:" << owner[faceI]
                    << " neighbour:" << ccFace.cells[1]
                    << exit(FatalError);
            }

            faceI++;
        }

    } // release the MeshCCMIO mesh.


    // Do some sanity checks
    forAll(region, faceI)
    {
        if
        (
            (neighbour[faceI] == -1 && region[faceI] == -1)
         || (neighbour[faceI] != -1 && region[faceI] != -1)
        )
        {
            FatalErrorIn("getFaceOrder") << "Problem"
                << " for face " << faceI << endl
                << "Face owner cell: " << owner[faceI]
                << ", face neighbour cell: " << neighbour[faceI]
                << " face region:" << region[faceI] << endl
                << abort(FatalError);
        }
    }


    // We now have extracted all info from CCMIO:
    // - coordinates (points)
    // - face to point addressing (faces)
    // - face to cell addressing (owner, neighbour, region)
    // - cell based data (cellId)
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

            if (own >= 0)
            {
                if (cellMap[own] == -1)
                {
                    cellMap[own] = nCells++;
                }
                owner[faceI] = cellMap[own];
            }

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

        // Compact cellId
        labelList newCellId(nCells, -1);

        forAll(cellId, oldCellI)
        {
            if (cellMap[oldCellI] != -1)
            {
                newCellId[cellMap[oldCellI]] = cellId[oldCellI];
            }
        }
        cellId.transfer(newCellId);

        Info<< "Compacted cells from " << cellId.size() << " down to "
            << nCells << endl;
    }



    // Create cells (inverse of face-to-cell addressing)
    cellList cells(nCells);

    {
        // Create cell-faces addressing.
        List<DynamicList<label,1,2> > cellFaces(cells.size());

        // Presize cellFaces for quads.
        forAll(cellFaces, cellI)
        {
            cellFaces[cellI].setSize(4);
        }

        forAll(owner, faceI)
        {
            if (owner[faceI] >= 0)
            {
                cellFaces[owner[faceI]].append(faceI);
            }
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

        if (nbr >= cells.size() || own >= cells.size())
        {
            FatalErrorIn(args.executable())
                << "face:" << faceI
                << " nbr:" << nbr
                << " own:" << own
                << exit(FatalError);
        }

        if (nbr >= 0)
        {
            if (nbr < own)
            {
                owner[faceI] = neighbour[faceI];
                neighbour[faceI] = own;
                faces[faceI] = faces[faceI].reverseFace();
            }
        }
    }


    // Determine face order for upper-triangular ordering and internal/external
    // face ordering
    labelList patchSizes(patchNames.size(), 0);


    labelList oldToNew
    (
        getFaceOrder
        (
            cells,
            owner,
            neighbour,
            region,
            patchNames.size(),
            patchSizes
        )
    );

    Info<< "Size per patch:" << patchSizes << endl;

    // Reorder faces accordingly
    reorderFaces(oldToNew, cells, faces, owner, neighbour);

    // Construct fvMesh (without patches)
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        points,
        faces,
        cells
    );

    // Create patches. Use patch types to determine what Foam types to generate.
    List<polyPatch*> newPatches(patchNames.size());

    label meshFaceI = nInternalFaces;

    forAll(newPatches, patchI)
    {
        const word& patchName  = patchNames[patchI];
        const word& patchType = patchTypes[patchI];

        if (patchType == "wall")
        {
            newPatches[patchI] =
                new wallPolyPatch
                (
                    patchName,
                    patchSizes[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh()
                );
        }
        else if (patchType == "symmetry")
        {
            newPatches[patchI] =
                new symmetryPolyPatch
                (
                    patchName,
                    patchSizes[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh()
                );
        }
        else if (patchType == "empty")
        {
            // Note: not ccm name, introduced by us above.
            newPatches[patchI] =
                new emptyPolyPatch
                (
                    patchName,
                    patchSizes[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh()
                );
        }
        else
        {
            // All other ccm types become straight polyPatch:
            // 'inlet', 'outlet', 'pressured'.
            newPatches[patchI] =
                new polyPatch
                (
                    patchName,
                    patchSizes[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh()
                );
        }

        meshFaceI += patchSizes[patchI];
    }

    mesh.addFvPatches(newPatches);

    //merge couples


    Info<< "Writing polyMesh to " << mesh.objectRegistry::objectPath()
        << "..." << nl << endl;

    mesh.write();


    // Construct field with calculated bc to hold Star cell Id.
    volScalarField cellIdField
    (
        IOobject
        (
            "cellId",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellId", dimless, 0.0)
    );

    forAll(cellId, cellI)
    {
        cellIdField[cellI] = cellId[cellI];
    }

    Info<< "Writing cellIds as volScalarField to " << cellIdField.objectPath()
        << "..." << nl << endl;
    cellIdField.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
