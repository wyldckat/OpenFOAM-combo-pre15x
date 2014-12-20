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

#include "edgeCollapser.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "ListSearch.H"
#include "polyAddFace.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::edgeCollapser::findIndex
(
    const labelList& elems,
    const label nElems,
    const label val
)
{
    for (label i = 0; i < nElems; i++)
    {
        if (elems[i] == val)
        {
            return i;
        }
    }
    return -1;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Changes region of connected set of points
Foam::label Foam::edgeCollapser::changePointRegion
(
    const label pointI,
    const label oldRegion,
    const label newRegion
)
{
    label nChanged = 0;

    if (pointRegion_[pointI] == oldRegion)
    {
        pointRegion_[pointI] = newRegion;
        nChanged++;

        // Step to neighbouring point across edges with same region number

        const labelList& pEdges = mesh_.pointEdges()[pointI];

        forAll(pEdges, i)
        {
            label otherPointI = mesh_.edges()[pEdges[i]].otherVertex(pointI);

            nChanged += changePointRegion(otherPointI, oldRegion, newRegion);
        }
    }
    return nChanged;
}


bool Foam::edgeCollapser::pointRemoved(const label pointI) const
{
    label region = pointRegion_[pointI];

    if (region == -1 || pointRegionMaster_[region] == pointI)
    {
        return false;
    }
    else
    {
        return true;
    }
}


void Foam::edgeCollapser::filterFace(face& f) const
{
    label newFp = 0;

    forAll(f, fp)
    {
        label pointI = f[fp];

        label region = pointRegion_[pointI];

        if (region == -1)
        {
            f[newFp++] = pointI;
        }
        else
        {
            label master = pointRegionMaster_[region];

            if (findIndex(f, newFp, master) == -1)
            {
                f[newFp++] = master;
            }
        }
    }
    f.setSize(newFp);
}


// Debugging.
void Foam::edgeCollapser::printRegions() const
{
    forAll(pointRegionMaster_, regionI)
    {
        label master = pointRegionMaster_[regionI];

        if (master != -1)
        {
            Info<< "Region:" << regionI << nl
                << "    master:" << master
                << ' ' << mesh_.points()[master] << nl;

            forAll(pointRegion_, pointI)
            {
                if (pointRegion_[pointI] == regionI && pointI != master)
                {
                    Info<< "    slave:" << pointI
                        << ' ' <<  mesh_.points()[pointI] << nl;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::edgeCollapser::edgeCollapser(const polyMesh& mesh)
:
    mesh_(mesh),
    pointRegion_(mesh.nPoints(), -1),
    pointRegionMaster_(mesh.nPoints(), -1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::edgeCollapser::unaffectedEdge(const label edgeI) const
{
    const edge& e = mesh_.edges()[edgeI];

    return (pointRegion_[e[0]] == -1) && (pointRegion_[e[1]] == -1);
}


bool Foam::edgeCollapser::collapseEdge(const label edgeI, const label master)
{
    const edge& e = mesh_.edges()[edgeI];

    label pointRegion0 = pointRegion_[e[0]];
    label pointRegion1 = pointRegion_[e[1]];

    if (pointRegion0 == -1)
    {
        if (pointRegion1 == -1)
        {
            // Both endpoints not collapsed. Create new region.
            label freeRegion =
                findIndex
                (
                    pointRegionMaster_,
                    pointRegionMaster_.size(),
                    -1
                );

            pointRegion_[e[0]] = freeRegion;
            pointRegion_[e[1]] = freeRegion;

            pointRegionMaster_[freeRegion] = master;
        }
        else
        {
            // e[1] is part of collapse network, e[0] not. Add e0 to e1 region.
            pointRegion_[e[0]] = pointRegion1;

            pointRegionMaster_[pointRegion1] = master;
        }
    }
    else
    {
        if (pointRegion1 == -1)
        {
            // e[0] is part of collapse network. Add e1 to e0 region
            pointRegion_[e[1]] = pointRegion0;

            pointRegionMaster_[pointRegion0] = master;
        }
        else if (pointRegion0 != pointRegion1)
        {
            // Both part of collapse network. Merge the two regions.

            // Use the smaller region number for the whole network.
            label minRegion = min(pointRegion0, pointRegion1);
            label maxRegion = max(pointRegion0, pointRegion1);
    
            // Use minRegion as region for combined net, free maxRegion.
            pointRegionMaster_[minRegion] = master;
            pointRegionMaster_[maxRegion] = -1;

            if (minRegion != pointRegion0)
            {
                changePointRegion(e[0], pointRegion0, minRegion);
            }
            if (minRegion != pointRegion1)
            {
                changePointRegion(e[1], pointRegion1, minRegion);
            }
        }
    }

    return true;
}


void Foam::edgeCollapser::setRefinement(polyTopoChange& meshMod) const
{
    const cellList& cells = mesh_.cells();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const labelListList& pointFaces = mesh_.pointFaces();


    // Print regions:
    //printRegions()


    // Current faces (is also collapseStatus: f.size() < 3)
    faceList newFaces(mesh_.faces());

    // Update face collapse from edge collapses
    forAll(newFaces, faceI)
    {
        filterFace(newFaces[faceI]);
    }

    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];

        label nFaces = cFaces.size();

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (newFaces[faceI].size() < 3)
            {
                --nFaces;

                if (nFaces < 4)
                {
                    Info<< "Cell:" << cellI
                        << " uses faces:" << cFaces
                        << " of which too many are marked for removal:" << endl
                        << "   ";
                    forAll(cFaces, j)
                    {
                        if (newFaces[cFaces[j]].size() < 3)
                        {
                            Info<< ' '<< cFaces[j];
                        }
                    }
                    Info<< endl;
                    FatalErrorIn("edgeCollapser::setRefinement(polyTopoChange&")
                        << "cell " << cellI << " would get less than 4 faces"
                        << abort(FatalError);
                }
            }
        }
    }


    // Keep track of faces that have been done already.
    boolList doneFace(mesh_.nFaces(), false);

    // Remove faces
    forAll(newFaces, faceI)
    {
        if (newFaces[faceI].size() < 3)
        {
            meshMod.setAction(polyRemoveFace(faceI));

            // Mark face as been done.
            doneFace[faceI] = true;
        }
    }


    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const faceZoneMesh& faceZones = mesh_.faceZones();
      

    // Renumber faces that use points
    forAll(pointRegion_, pointI)
    {
        if (pointRemoved(pointI))
        {
            const labelList& changedFaces = pointFaces[pointI];

            forAll(changedFaces, changedFaceI)
            {
                label faceI = changedFaces[changedFaceI];

                if (!doneFace[faceI])
                {
                    doneFace[faceI] = true;

                    // Get current zone info
                    label zoneID = faceZones.whichZone(faceI);

                    bool zoneFlip = false;

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = faceZones[zoneID];

                        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                    }

                    // Get current connectivity
                    label own = faceOwner[faceI];
                    label nei = -1;
                    label patchID = -1;

                    if (mesh_.isInternalFace(faceI))
                    {
                        nei = faceNeighbour[faceI];
                    }
                    else
                    {
                        patchID = boundaryMesh.whichPatch(faceI);
                    }

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            newFaces[faceI],            // face
                            faceI,                      // faceI to change
                            own,                        // owner
                            nei,                        // neighbour
                            false,                      // flipFaceFlux
                            patchID,                    // patch
                            false,                      // removeFromZone
                            zoneID,                     // zoneID
                            zoneFlip                    // zoneFlip
                        )
                    );
                }
            }
        }
    }
}


// ************************************************************************* //
