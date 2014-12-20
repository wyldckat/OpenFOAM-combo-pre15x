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

#include "directPolyTopoChange.H"
#include "polyMesh.H"
#include "Time.H"
#include "morphMesh.H"
#include "SortableList.H"
#include "PackedList.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directPolyTopoChange, 0);
}



const Foam::point Foam::directPolyTopoChange::greatPoint
(
    GREAT,
    GREAT,
    GREAT
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Renumber
void Foam::directPolyTopoChange::renumber
(
    const labelList& map,
    DynamicList<label>& elems
)
{
    forAll(elems, elemI)
    {
        if (elems[elemI] != -1)
        {
            elems[elemI] = map[elems[elemI]];
        }
    }
}


void Foam::directPolyTopoChange::remove
(
    const label elemToRemove,
    labelList& elems
)
{
    label newElemI = 0;

    forAll(elems, elemI)
    {
        if (elems[elemI] != elemToRemove)
        {
            elems[newElemI++] = elems[elemI];
        }
    }
    elems.setSize(newElemI);
}
    

// Renumber and remove -1 elements.
void Foam::directPolyTopoChange::renumberCompact
(
    const labelList& map,
    labelList& elems
)
{
    label newElemI = 0;

    forAll(elems, elemI)
    {
        label newVal = map[elems[elemI]];

        if (newVal != -1)
        {
            elems[newElemI++] = newVal;
        }
    }
    elems.setSize(newElemI);
}


bool Foam::directPolyTopoChange::pointRemoved(const point& pt)
{
    return
        pt.x() > 0.5*greatPoint.x()
     && pt.y() > 0.5*greatPoint.y()
     && pt.z() > 0.5*greatPoint.z();
}


// Determine order for faces:
// - upper-triangular order for internal faces
// - external faces after internal faces
Foam::labelList Foam::directPolyTopoChange::getFaceOrder(const label nPatches)
 const
{
    labelList oldToNew(faceOwner_.size(), -1);

    // First unassigned face
    label newFaceI = 0;

    forAll(cells_, cellI)
    {
        const labelList& cFaces = cells_[cellI];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            label nbrCellI = faceNeighbour_[faceI];

            if (nbrCellI != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCellI == cellI)
                {
                    nbrCellI = faceOwner_[faceI];
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
    for (label patchI = 0; patchI < nPatches; patchI++)
    {
        forAll(regions_, faceI)
        {
            if (regions_[faceI] == patchI)
            {
                oldToNew[faceI] = newFaceI++;
            }
        }
    }


    // Check done all faces.
    forAll(oldToNew, faceI)
    {
        if (oldToNew[faceI] == -1)
        {
            FatalErrorIn
            (
                "directPolyTopoChange::getFaceOrder"
                "(const labelList&, const labelList&, const label) const"
            )   << "Did not determine new position"
                << " for face " << faceI
                << abort(FatalError);
        }
    }

    return oldToNew;
}


void Foam::directPolyTopoChange::compact()
{
    points_.shrink();
    faces_.shrink();
    regions_.shrink();
    faceOwner_.shrink();
    faceNeighbour_.shrink();
    cells_.shrink();

    pointMap_.shrink();
    faceMap_.shrink();
    cellMap_.shrink();

    if (debug)
    {
        Info<< "Before compacting:" << nl
            << "    points:" << points_.size() << nl
            << "    faces :" << faces_.size() << nl
            << "    cells :" << cells_.size() << nl
            << endl;
    }

    // Remove unused faces
    {
        PackedList<1> usedFace(faces_.size(), 0);

        forAll(cells_, cellI)
        {
            const cell& cFaces = cells_[cellI];

            forAll(cFaces, i)
            {
                label faceI = cFaces[i];

                if (faces_[faceI].size() == 0)
                {
                    FatalErrorIn("directPolyTopoChange::compact()")
                        << "Cell " << cellI << " faces:" << cFaces
                        << " uses removed face " << faceI
                        << abort(FatalError);
                }

                usedFace.set(faceI, 1);
            }
        }

        forAll(usedFace, faceI)
        {
            if (!usedFace.get(faceI) && faces_[faceI].size() != 0)
            {
                removeFace(faceI);
            }
        }
    }

    // Remove unused points
    {
        PackedList<1> usedPoint(points_.size(), 0);

        forAll(faces_, faceI)
        {
            const face& f = faces_[faceI];

            forAll(f, fp)
            {
                label pointI = f[fp];

                if (pointRemoved(points_[pointI]))
                {
                    FatalErrorIn("directPolyTopoChange::compact()")
                        << "Face " << faceI << " vertices:" << f
                        << " uses removed point " << pointI
                        << abort(FatalError);
                }

                usedPoint.set(f[fp], 1);
            }
        }

        forAll(usedPoint, pointI)
        {
            if (!usedPoint.get(pointI) && !pointRemoved(points_[pointI]))
            {
                removePoint(pointI);
            }
        }
    }

    if (debug)
    {
        Info<< "After removing unused faces and points:" << nl
            << "    points:" << points_.size() << nl
            << "    faces :" << faces_.size() << nl
            << "    cells :" << cells_.size() << nl
            << endl;
    }

    // Compact points
    {
        labelList localPointMap(points_.size(), -1);
        label newPointI = 0;

        forAll(points_, pointI)
        {
            if (!pointRemoved(points_[pointI]))
            {
                localPointMap[pointI] = newPointI++;
            }
        }

        if (newPointI != points_.size())
        {
            reorder(localPointMap, points_);
            points_.setSize(newPointI);
            points_.shrink();

            // Use map to relabel face vertices
            forAll(faces_, faceI)
            {
                face& f = faces_[faceI];

                labelList oldF(f);

                renumberCompact(localPointMap, f);

                if (f.size() != 0 && f.size() < 3)
                {
                    FatalErrorIn("directPolyTopoChange::compact()")
                        << "Created illegal face " << f
                        << " from face " << oldF
                        << " at position:" << faceI
                        << " when filtering removed points"
                        << abort(FatalError);
                }
            }

            // Update overall pointMap
            renumber(localPointMap, pointMap_);
        }
    }


    // Compact faces.
    {
        labelList localFaceMap(faces_.size(), -1);
        label newFaceI = 0;

        forAll(faces_, faceI)
        {
            const face& f = faces_[faceI];

            if (f.size() != 0)
            {
                localFaceMap[faceI] = newFaceI++;
            }
        }

        if (newFaceI != faces_.size())
        {
            reorder(localFaceMap, faces_);
            faces_.setSize(newFaceI);
            faces_.shrink();

            reorder(localFaceMap, regions_);
            regions_.setSize(newFaceI);
            regions_.shrink();

            reorder(localFaceMap, faceOwner_);
            faceOwner_.setSize(newFaceI);
            faceOwner_.shrink();

            reorder(localFaceMap, faceNeighbour_);
            faceNeighbour_.setSize(newFaceI);
            faceNeighbour_.shrink();

            forAll(cells_, cellI)
            {
                cell& c = cells_[cellI];

                labelList oldC(c);

                renumberCompact(localFaceMap, c);

                if (c.size() != 0 && c.size() < 4)
                {
                    FatalErrorIn("directPolyTopoChange::compact()")
                        << "Created illegal cell " << c
                        << " from cell " << oldC
                        << " at position:" << cellI
                        << " when filtering removed faces"
                        << abort(FatalError);
                }
            }

            // Update overall faceMap.
            renumber(localFaceMap, faceMap_);
        }
    }

    // Compact cells.
    {
        labelList localCellMap(cells_.size(), -1);
        label newCellI = 0;

        forAll(cells_, cellI)
        {
            if (cells_[cellI].size() != 0)
            {
                localCellMap[cellI] = newCellI++;
            }
        }

        if (newCellI != cells_.size())
        {
            reorder(localCellMap, cells_);
            cells_.setSize(newCellI);
            cells_.shrink();

            renumber(localCellMap, faceOwner_);
            renumber(localCellMap, faceNeighbour_);

            // Update overall cellMap.
            renumber(localCellMap, cellMap_);
        }
    }

    // Reorder faces into upper-triangular and patch ordering
    {
        // Count regions
        label maxRegion = -1;
        forAll(regions_, faceI)
        {
            maxRegion = max(maxRegion, regions_[faceI]);
        }
        label nPatches = maxRegion + 1;

        // Do upper triangular order.
        labelList localFaceMap(getFaceOrder(nPatches));

        reorder(localFaceMap, faces_);
        reorder(localFaceMap, regions_);
        reorder(localFaceMap, faceOwner_);
        reorder(localFaceMap, faceNeighbour_);

        forAll(cells_, cellI)
        {
            renumberCompact(localFaceMap, cells_[cellI]);
        }

        // Update overall face map.
        renumber(localFaceMap, faceMap_);
    }

    if (debug)
    {
        Info<< "After removing removed points,faces,cell:" << nl
            << "    points:" << points_.size() << nl
            << "    faces :" << faces_.size() << nl
            << "    cells :" << cells_.size() << nl
            << endl;
    }
}


void Foam::directPolyTopoChange::calcPatchSizes
(
    label& nInternalFaces,
    labelList& patchSizes
) const
{
    label maxRegion = -1;
    forAll(regions_, faceI)
    {
        maxRegion = max(maxRegion, regions_[faceI]);
    }
    label nPatches = maxRegion + 1;

    nInternalFaces = 0;

    patchSizes.setSize(nPatches);
    patchSizes = 0;

    forAll(regions_, faceI)
    {
        if (regions_[faceI] != -1)
        {
            patchSizes[regions_[faceI]]++;
        }
        else
        {
            nInternalFaces++;
        } 
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::directPolyTopoChange::directPolyTopoChange(const bool strict)
:
    points_(),
    faces_(),
    regions_(),
    faceOwner_(),
    faceNeighbour_(),
    cells_(),
    pointMap_(),
    faceMap_(),
    cellMap_(),
    strict_(strict)
{}


// Construct from components
Foam::directPolyTopoChange::directPolyTopoChange
(
    const polyMesh& mesh,
    const bool strict
)
:
    points_(mesh.nPoints()),
    faces_(mesh.nFaces()),
    regions_(mesh.nFaces()),
    faceOwner_(mesh.nFaces()),
    faceNeighbour_(mesh.nFaces()),
    cells_(mesh.nCells()),
    pointMap_(mesh.nPoints()),
    faceMap_(mesh.nFaces()),
    cellMap_(mesh.nCells()),
    strict_(strict)
{
    // Do like addMesh but copy whole cells in one go.

    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    // Add points in mesh order
    forAll(points, pointI)
    {
        addPoint(points[pointI]);
    }

    // Add faces in mesh order
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        faces_.append(faces[faceI]);
        regions_.append(-1);
        faceOwner_.append(faceOwner[faceI]);
        faceNeighbour_.append(faceNeighbour[faceI]);
        faceMap_.append(faceI);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.start() != faces_.size())
        {
            FatalErrorIn
            (
                "directPolyTopoChange::directPolyTopoChange"
                "(const polyMesh& mesh, const bool strict)"
            )   << "Problem : "
                << "Patch " << pp.name() << " starts at " << pp.start() << endl
                << "Current face counter at " << faces_.size() << endl
                << "Are patches in incremental order?"
                << abort(FatalError);
        }
        forAll(pp, patchFaceI)
        {
            label faceI = pp.start() + patchFaceI;

            faces_.append(faces[faceI]);
            regions_.append(patchI);
            faceOwner_.append(faceOwner[faceI]);
            faceNeighbour_.append(-1);
            faceMap_.append(faceI);
        }
    }

    // Add cells
    forAll(cells, cellI)
    {
        cells_.append(cells[cellI]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directPolyTopoChange::addMesh(const polyMesh& mesh)
{
    // points
    const pointField& points = mesh.points();

    forAll(points, pointI)
    {
        addPoint(points[pointI]);
    }


    // cells
    const cellList& cells = mesh.cells();

    forAll(cells, cellI)
    {
        addCell();
    }


    // faces
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        addFace(faces[faceI], faceOwner[faceI], faceNeighbour[faceI], -1);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll(pp, patchFaceI)
        {
            addFace
            (
                pp[patchFaceI],
                faceOwner[pp.start() + patchFaceI],
                -1,
                patchI
            );
        }
    }
}


Foam::label Foam::directPolyTopoChange::addPoint(const point& pt)
{
    label pointI = points_.size();

    points_.append(pt);
    pointMap_.append(pointI);

    return pointI;
}


void Foam::directPolyTopoChange::modifyPoint
(
    const label pointI,
    const point& pt
)
{
    if (pointI < 0 || pointI >= points_.size())
    {
        FatalErrorIn
        (
            "directPolyTopoChange::modifyPoint(const label, const point&)"
        )   << "illegal point label " << pointI << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }
    if (pointRemoved(points_[pointI]) || pointMap_[pointI] == -1)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::modifyPoint(const label, const point&)"
        )   << "point " << pointI << " already marked for removal"
            << abort(FatalError);
    }
    points_[pointI] = pt;
}


void Foam::directPolyTopoChange::removePoint(const label pointI)
{
    if (pointI < 0 || pointI >= points_.size())
    {
        FatalErrorIn("directPolyTopoChange::removePoint(const label)")
            << "illegal point label " << pointI << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (pointRemoved(points_[pointI]) || pointMap_[pointI] == -1)
    )
    {
        FatalErrorIn("directPolyTopoChange::removePoint(const label)")
            << "point " << pointI << " already marked for removal" << nl
            << "Point:" << points_[pointI] << " pointMap:" << pointMap_[pointI]
            << abort(FatalError);
    }

    points_[pointI] = greatPoint;
    pointMap_[pointI] = -1;
}


Foam::label Foam::directPolyTopoChange::addFace
(
    const face& f,
    const label own,
    const label nei,
    const label patchI
)
{
    if (patchI != -1 && nei != -1)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::addFace(const face&, const label"
            ", const label, const label)"
        )   << "Cannot both have valid patchI and neighbour" << nl
            << "f:" << f << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << abort(FatalError);
    }

    if (nei != -1 && nei <= own)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::addFace(const face&, const label"
            ", const label, const label)"
        )   << "Owner cell label should be less than neighbour cell label"
            << nl
            << "f:" << f << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << abort(FatalError);
    }

    label faceI = faces_.size();


    faces_.append(f);
    regions_.append(patchI);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);
    faceMap_.append(faceI);

    if (own >= cells_.size())
    {
        if (strict_)
        {
            FatalErrorIn
            (
                "directPolyTopoChange::addFace(const face&, const label"
                ", const label, const label)"
            )   << "Illegal owner cell label" << nl
                << "f:" << f << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << abort(FatalError);
        }

        for (label cellI = cells_.size(); cellI <= own; cellI++)
        {
            // Add dummy cell to make sure cellMap is set.
            addCell();
        }
    }
    cell& cFaces = cells_[own];

    label i = cFaces.size();
    cFaces.setSize(i+1);
    cFaces[i] = faceI;

    if (nei != -1)
    {
        if (nei >= cells_.size())
        {
            if (strict_)
            {
                FatalErrorIn
                (
                    "directPolyTopoChange::addFace(const face&, const label"
                    ", const label, const label)"
                )   << "Illegal neighbour cell label" << nl
                    << "f:" << f << " own:" << own << " nei:" << nei
                    << " patchI:" << patchI << abort(FatalError);
            }

            for (label cellI = cells_.size(); cellI <= nei; cellI++)
            {
                // Add dummy cell to make sure cellMap is set.
                addCell();
            }
        }
        cell& cFaces = cells_[nei];

        label i = cFaces.size();
        cFaces.setSize(i+1);
        cFaces[i] = faceI;
    }

    return faceI;
}    


void Foam::directPolyTopoChange::modifyFace
(
    const face& f,
    const label faceI,
    const label own,
    const label nei,
    const label patchI
)
{
    if (patchI != -1 && nei != -1)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::addFace(const face&, const label"
            ", const label, const label)"
        )   << "Cannot both have valid patchI and neighbour" << nl
            << "f:" << f << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << abort(FatalError);
    }

    if (nei != -1 && nei <= own)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::addFace(const face&, const label"
            ", const label, const label)"
        )   << "Owner cell label should be less than neighbour cell label"
            << nl
            << "f:" << f << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << abort(FatalError);
    }

    if (faceOwner_[faceI] != -1 && own != faceOwner_[faceI])
    {
        // Remove faceI from faceOwner cell
        remove(faceI, cells_[faceOwner_[faceI]]);
    }
    if (faceNeighbour_[faceI] != -1 && nei != faceNeighbour_[faceI])
    {
        // Remove faceI from faceNeighbour cell
        remove(faceI, cells_[faceNeighbour_[faceI]]);
    }

    faces_[faceI] = f;
    regions_[faceI] = patchI;
    faceMap_[faceI] = faceI;

    if (own != faceOwner_[faceI])
    {
        cell& cFaces = cells_[own];

        label i = cFaces.size();
        cFaces.setSize(i+1);
        cFaces[i] = faceI;
        faceOwner_[faceI] = own;
    }

    if (nei != -1 && nei != faceNeighbour_[faceI])
    {
        cell& cFaces = cells_[nei];

        label i = cFaces.size();
        cFaces.setSize(i+1);
        cFaces[i] = faceI;
        faceNeighbour_[faceI] = nei;
    }
}    


void Foam::directPolyTopoChange::removeFace(const label faceI)
{
    if (faceI < 0 || faceI >= faces_.size())
    {
        FatalErrorIn("directPolyTopoChange::removeFace(const label)")
            << "illegal face label " << faceI << endl
            << "Valid face labels are 0 .. " << faces_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (faces_[faceI].size() == 0 || faceMap_[faceI] == -1)
    )
    {
        FatalErrorIn("directPolyTopoChange::removeFace(const label)")
            << "face " << faceI
            << " already marked for removal"
            << abort(FatalError);
    }


    // Remove faceI from cells using it.
    if (faceOwner_[faceI] != -1)
    {
        remove(faceI, cells_[faceOwner_[faceI]]);
    }
    if (faceNeighbour_[faceI] != -1)
    {
        remove(faceI, cells_[faceNeighbour_[faceI]]);
    }

    faces_[faceI].setSize(0);
    faceOwner_[faceI] = -1;
    faceNeighbour_[faceI] = -1;
    faceMap_[faceI] = -1;
}


Foam::label Foam::directPolyTopoChange::addCell()
{
    // Append dummy cell
    label cellI = cells_.size();

    cells_.append(cell());
    cellMap_.append(cellI);

    return cellI;
}


void Foam::directPolyTopoChange::removeCell(const label cellI)
{
    if (cellI < 0 || cellI >= cells_.size())
    {
        FatalErrorIn("directPolyTopoChange::removeCell(const label)")
            << "illegal cell label " << cellI << endl
            << "Valid cell labels are 0 .. " << cells_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (cells_[cellI].size() == 0 || cellMap_[cellI] == -1)
    )
    {
        FatalErrorIn("directPolyTopoChange::removeCell(const label)")
            << "cell " << cellI
            << " already marked for removal"
            << abort(FatalError);
    }

    const cell& cFaces = cells_[cellI];

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];

        if (faceOwner_[faceI] == cellI)
        {
            faceOwner_[faceI] = -1;
        }
        else if (faceNeighbour_[faceI] == cellI)
        {
            faceNeighbour_[faceI] = -1;
        }
    }

    cells_[cellI].setSize(0);
    cellMap_[cellI] = -1;
}


Foam::autoPtr<Foam::morphMesh> Foam::directPolyTopoChange::mesh
(
    const IOobject& io
)
{
    if (debug)
    {
        Info<< "directPolyTopoChange::mesh without patches" << endl;
    }

    compact();

    // Transfer points to pointField
    pointField newPoints;
    newPoints.transfer(points_);
    points_.clear();

    // Construct mesh
    autoPtr<morphMesh> meshPtr
    (
        new morphMesh
        (
            io,
            newPoints,
            faces_,
            cells_
        )
    );

    // Clear all we don't use anymore
    faces_.faceList::clear();
    faces_.clear();

    faceOwner_.labelList::clear();
    faceOwner_.clear();
    faceNeighbour_.labelList::clear();
    faceNeighbour_.clear();

    cells_.cellList::clear();
    cells_.clear();

    // Leave regions_ intact.

    if (debug)
    {
        Info<< "directPolyTopoChange::mesh : created morphMesh" << endl;
    }

    return meshPtr;
}


// Sort faces on coupled patches to be consistent.
// From polyMeshMorph.C and couplePatches utility
void Foam::directPolyTopoChange::reorderCoupledPatches(morphMesh& mesh)
{
    if (debug)
    {
        Info<< "directPolyTopoChange::reorderCoupledPatches" << endl;
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    bool hasCoupled = false;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            hasCoupled = true;

            break;
        }
    }

    reduce(hasCoupled, orOp<bool>());

    if (hasCoupled)
    {
        if (debug)
        {
            Info<< "Mesh has coupled patches ..." << nl << endl;

            Info<< "Testing for correct face ordering ..." << endl;
        }

        bool wrongOrder = mesh.boundaryMesh().checkFaceOrder(debug);

        if (wrongOrder)
        {
            // Dummy topo changes container
            polyTopoChange meshMod(mesh);

            // Do all changes
            if (debug)
            {
                Info<< "Doing dummy mesh morph to correct face ordering ..."
                    << endl;
            }

            mesh.setMorphTimeIndex(mesh.time().timeIndex());

            mesh.updateTopology(meshMod);

            // Move mesh (since morphing does not do this)
            mesh.movePoints(mesh.morphMap().preMotionPoints());

            // Update maps.
            renumber(mesh.morphMap().reverseCellMap(), cellMap_);
            renumber(mesh.morphMap().reverseFaceMap(), faceMap_);
            renumber(mesh.morphMap().reversePointMap(), pointMap_);
        }
        else
        {
            if (debug)
            {
                Info<< "Coupled patch face ordering ok. Nothing changed ..."
                    << nl << endl;
            }
        }
    }
}


Foam::autoPtr<Foam::morphMesh> Foam::directPolyTopoChange::mesh
(
    const wordList& oldPatchNames,
    const wordList& oldPatchTypes,
    const IOobject& io
)
{
    if (debug)
    {
        Info<< "directPolyTopoChange::mesh with patches" << endl;
    }

    if (oldPatchNames.size() != oldPatchTypes.size())
    {
        FatalErrorIn
        (
            "directPolyTopoChange::mesh(const wordList&, const wordList&"
            ", const IOobject& io)"
        )   << "Sizes of oldPatchNames and oldPatchTypes not equal" << endl
            << "oldPatchNames:" << oldPatchNames
            << "  oldPatchTypes:" << oldPatchTypes << abort(FatalError);
    }

    // Get mesh without patches. (note: clears out cells_, points_, faces_)
    autoPtr<morphMesh> meshPtr = mesh(io);

    // Get patch sizes
    labelList patchSizes;
    label nInternalFaces = -1;
    calcPatchSizes(nInternalFaces, patchSizes);

    // Can clear out regions_ now.
    regions_.labelList::clear();
    regions_.clear();

    if (debug)
    {
        Info<< "nInternalFaces:" << nInternalFaces << endl;
        Info<< "patchSizes:" << patchSizes << endl;
    }

    List<polyPatch*> newPatches(max(oldPatchNames.size(), patchSizes.size()));

    label endOfLastPatch = nInternalFaces;

    // Copy all old patches (even if size 0)
    forAll(oldPatchNames, patchI)
    {
        // See if any faces in patch
        label sz = 0;

        if (patchI < patchSizes.size())
        {
            sz = patchSizes[patchI];
        }

        if (debug)
        {
            Info<< "adding patch " << patchI << " size:" << sz
                << " start:" << endOfLastPatch << endl;
        }

        newPatches[patchI] =
        (
            polyPatch::New
            (
                oldPatchTypes[patchI],
                oldPatchNames[patchI],
                sz,
                endOfLastPatch,
                patchI,
                meshPtr().boundaryMesh()
            ).ptr()
        );
        endOfLastPatch += sz;
    }

    // Add any extra patches
    for
    (
        label patchI = oldPatchNames.size();
        patchI < patchSizes.size();
        patchI++
    )
    {
        if (debug)
        {
            Info<< "adding extra patch " << patchI
                << " size:" << patchSizes[patchI]
                << " start:" << endOfLastPatch << endl;
        }

        newPatches[patchI] =
        (
            polyPatch::New
            (
                polyPatch::typeName,
                "patch" + name(patchI),
                patchSizes[patchI],
                endOfLastPatch,
                patchI,
                meshPtr().boundaryMesh()
            ).ptr()
        );
        endOfLastPatch += patchSizes[patchI];
    }

    meshPtr().addPatches(newPatches);

    reorderCoupledPatches(meshPtr());

    if (debug)
    {
        Info<< "directPolyTopoChange::mesh : done whole mesh" << endl;
    }


    return meshPtr;    
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
