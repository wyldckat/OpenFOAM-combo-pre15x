/*---------------------------------------------------------------------------* \
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
    Change topology of a polyMesh and create mesh-to-mesh mapping information

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "Time.H"
#include "faceZone.H"
#include "polyTopoChange.H"
#include "pointField.H"
#include "polyAddPoint.H"
#include "polyModifyPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "Map.H"
#include "DynamicList.H"
#include "primitiveMesh.H"
#include "mapPolyMesh.H"
#include "objectMap.H"
#include "parallelInfo.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::face Foam::polyMesh::rotateFace(const face& f, const label nPos)
{
    face newF(f.size());

    forAll(f, fp)
    {
        label fp1 = (fp + nPos) % f.size();

        if (fp1 < 0)
        {
            fp1 += f.size();
        }

        newF[fp1] = f[fp];
    }

    return newF;
}


bool Foam::polyMesh::reorderCoupledPatches
(
    const polyTopoChange& ref,
    mapPolyMesh& morphMap
)
{
    // Mapping for faces (old to new). Extends over all mesh faces for
    // convenience (could be just the external faces)
    labelList faceMap(faces_.size());
    forAll(faceMap, faceI)
    {
        faceMap[faceI] = faceI;
    }

    // Rotation on new faces.
    labelList rotation(faces_.size(), 0);


    // Prepare patches for face reordering
    forAll (boundary_, patchI)
    {
        boundary_[patchI].initOrder(ref, morphMap);
    }

    // Get patch face reordering
    bool anyChanged = false;

    forAll (boundary_, patchI)
    {
        // Get ordering for patch
        labelList patchFaceMap;
        labelList patchFaceRotation;

        bool changed =
            boundary_[patchI].order
            (
                ref,
                *morphMap_,
                patchFaceMap,
                patchFaceRotation
            );

        if (changed)
        {
            // Merge patch face reordering into mesh face reordering table
            label start = boundary_[patchI].start();

            forAll(patchFaceMap, patchFaceI)
            {
                faceMap[patchFaceI + start] = start + patchFaceMap[patchFaceI];
            }

            forAll(patchFaceRotation, patchFaceI)
            {
                rotation[patchFaceI + start] = patchFaceRotation[patchFaceI];
            }

            anyChanged = true;
        }
    }

    reduce(anyChanged, orOp<bool>());

    if (anyChanged)
    {
        // Update anything to do with face labels:
        //  - xxxFace and faceXXX addressing:
        //      - cell-face (cells)
        //      - face-face
        //      - edge-face
        //      - point-face
        //      - face-point (faces)
        //      - face-edge
        //      - face-cell (owner)
        //  - patches (done above)
        //  - faceZones
        //  - morphMap

        // Working copies
        faceList newFaces(faces_);
        cellList newCells(cells_);

        forAll(faceMap, faceI)
        {
            label newFaceI = faceMap[faceI];

            if (newFaceI != faceI)
            {
                // Face gets remapped.

                // face-point
                newFaces[newFaceI] =
                    rotateFace(faces_[faceI], rotation[newFaceI]);

                // cell-face
                label own = faceOwner()[faceI];

                // Index of old face in old cell.
                label index = findIndex(cells_[own], faceI);

                newCells[own][index] = newFaceI;
            }
            else if (rotation[newFaceI] != 0)
            {
                // Face label stays same but face gets rotated

                newFaces[newFaceI] =
                    rotateFace(faces_[faceI], rotation[newFaceI]);
            }   
        }

        // Clear all topological addressing. All geometry stays same.
        clearAddressing();
        clearFaceCells();

        // Set new faces and cells.
        faces_ = newFaces;
        cells_ = newCells;

        // Redo owner/neighbour calc.
        calcFaceCells();

        // Reset the boundary patches
        forAll (boundary_, patchI)
        {
            boundary_[patchI] = polyPatch
            (
                boundary_[patchI].name(),
                boundary_[patchI].size(),
                boundary_[patchI].start(),
                patchI,
                boundary_
            );
        }

        // Update mapPolyMesh
        morphMap.reorderPatchFaces(faceMap);

        // Upate faceZones
        forAll (faceZones_, fzI)
        {
            faceZones_[fzI].updateTopology(morphMap);
        }
    }

    return anyChanged;
}


const Foam::polyMeshMorphEngine& Foam::polyMesh::morphEngine() const
{
    if (!morphEnginePtr_)
    {
        morphEnginePtr_ = new polyMeshMorphEngine
        (
            IOobject
            (
                "meshModifiers",
                time().findInstance
                (
                    meshDir(),
                    "meshModifiers",
                    IOobject::READ_IF_PRESENT
                ),
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this
        );
    }

    return *morphEnginePtr_;
}


void Foam::polyMesh::morph
(
    const polyTopoChange& ref
)
{
    // It is not allowed to call a topology change twice in the same
    // time-step because of the limitations of zero volume/area
    // object addition/removal.  

    // If the morph has already been called in this time-step,
    // the morph map pointer will be active (it is deleted in the
    // first morph check on time-step increment
    if (morphMap_)
    {
        FatalErrorIn
        (
            "void polyMesh::morph(const polyTopoChange& ref)"
        )   << "Requested second topological change within a "
            << "single time-step." << nl << "This is not allowed."
            << abort(FatalError);
    }

    // Create a new list of points
    // Note.  Modified points are only influenced in the mesh motion stage
    if (debug || morphDebug)
    {
        Pout<< "Foam::polyMesh::morph" << nl
            << '(' << nl
            << "    const polyTopoChange& ref" << nl
            << ") : started executing topological change." << nl << endl;

        // Check topological change request for consistency and report
        // the statistics of the refinement request
        if (ref.check())
        {
            FatalErrorIn
            (
                "void polyMesh::morph(const polyTopoChange& ref)"
            )   << "Inconsistent topological change request."
                << abort(FatalError);
        }
    }


    // See if we need to do any coupled patch face ordering afterwards.
    bool hasCoupled = false;

    forAll (boundary_, patchi)
    {
        hasCoupled |= boundary_[patchi].coupled();
    }
    reduce(hasCoupled, orOp<bool>());


    // Grab old live mesh sizes
    const label nOldPoints = nPoints();
    const label nOldFaces = nFaces();
    const label nOldCells = nCells();

    // Keep the old patch start labels
    labelList oldPatchStarts(boundary_.size());
    forAll (boundary_, patchi)
    {
        oldPatchStarts[patchi] = boundary_[patchi].start();
    }

    pointField newPointsZeroVol(points_.size() + ref.pointBalance());
    pointField newPointsMotion(points_.size() + ref.pointBalance());

    // renumberPoints contains the new point label for every old and added point
    labelList renumberPoints(points_.size() + ref.addedPoints().size(), -1);

    // pointMap contains the old point label or the master point label
    // for all new points
    labelList pointMap(points_.size() + ref.addedPoints().size(), -1);
    label nNewPoints = 0;

    const labelHashSet& removedPoints = ref.removedPoints();
    const labelHashSet& removedFaces = ref.removedFaces();
    const labelHashSet& removedCells = ref.removedCells();

    // Grab the untouched points.
    forAll (points_, pointI)
    {
        // Check if the point has been removed; if not add it to the list
        if (!removedPoints.found(pointI))
        {
            // Grab a point
            newPointsZeroVol[nNewPoints] = points_[pointI];
            newPointsMotion[nNewPoints] = points_[pointI];

            // Grab addressing
            renumberPoints[pointI] = nNewPoints;
            pointMap[nNewPoints] = pointI;

            nNewPoints++;
        }
    }

    label debugPointCounter = 0;

    if (morphDebug)
    {
        Pout<< "Added untouched points. Point count = "
            << nNewPoints << endl;

        debugPointCounter = nNewPoints;
    }

    // Change the modified points in two passes: first points
    // supporting the cells and then auxiliary points
    const DynamicList<polyModifyPoint>& mp = ref.modifiedPoints();

    forAll (mp, mpI)
    {
        if (mp[mpI].inCell())
        {
            // Grab a point
            newPointsZeroVol[nNewPoints] = points_[mp[mpI].pointID()];
            newPointsMotion[nNewPoints] = mp[mpI].newPoint();

            // Grab addressing
            renumberPoints[mp[mpI].pointID()] = nNewPoints;
            pointMap[nNewPoints] = mp[mpI].pointID();
            nNewPoints++;
        }
    }

    if (morphDebug)
    {
        Pout<< "Added live points: modified = "
            << nNewPoints - debugPointCounter;

        debugPointCounter = nNewPoints;
    }

    // Add the new points in two passes: first points supporting the cells
    // and then auxiliary points
    const DynamicList<polyAddPoint>& ap = ref.addedPoints();

    // Follow the location in the renumbering list
    const label np = points_.size();

    // Grab points supporting cells
    forAll (ap, apI)
    {
        if (ap[apI].inCell())
        {
            // Grab a point
            if (ap[apI].appended())
            {
                // Point appended directly; no motion
                newPointsZeroVol[nNewPoints] = ap[apI].newPoint();
            }
            else
            {
                // Point added on top of another point
                newPointsZeroVol[nNewPoints] = points_[ap[apI].masterPointID()];
            }

            newPointsMotion[nNewPoints] = ap[apI].newPoint();

            // Grab addressing
            renumberPoints[np + apI] = nNewPoints;
            pointMap[nNewPoints] = ap[apI].masterPointID();
            nNewPoints++;
        }
    }

    if (morphDebug)
    {
        Pout<< " added = " << nNewPoints - debugPointCounter
            << ".  Point count = "
            << nNewPoints << endl;

        debugPointCounter = nNewPoints;
    }

    // Grab auxiliary points
    forAll (mp, mpI)
    {
        if (!mp[mpI].inCell())
        {
            // Grab a point
            newPointsZeroVol[nNewPoints] = points_[mp[mpI].pointID()];
            newPointsMotion[nNewPoints] = mp[mpI].newPoint();

            // Grab addressing
            renumberPoints[mp[mpI].pointID()] = nNewPoints;
            pointMap[nNewPoints] = mp[mpI].pointID();
            nNewPoints++;
        }
    }

    if (morphDebug)
    {
        Pout<< "Added retired points: modified = "
            << nNewPoints - debugPointCounter;

        debugPointCounter = nNewPoints;
    }

    forAll (ap, apI)
    {
        if (!ap[apI].inCell())
        {
            // Grab a point
            if (ap[apI].appended())
            {
                // Point appended directly; no motion
                newPointsZeroVol[nNewPoints] = ap[apI].newPoint();
            }
            else
            {
                // Point added on top of another point
                newPointsZeroVol[nNewPoints] = points_[ap[apI].masterPointID()];
            }

            newPointsMotion[nNewPoints] = ap[apI].newPoint();

            // Grab addressing
            renumberPoints[np + apI] = nNewPoints;
            pointMap[nNewPoints] = ap[apI].masterPointID();
            nNewPoints++;
        }
    }

    if (morphDebug)
    {
        Pout<< " added = " << nNewPoints - debugPointCounter
            << ".  Point count = "
            << nNewPoints << endl;

        debugPointCounter = nNewPoints;
    }

    // Reset the size of point map
    pointMap.setSize(nNewPoints);

    if (morphDebug)
    {
        Pout<< "Added all points.  Final point count = "
            << nNewPoints << nl << endl;
    }

    // Create a new list of faces
    // Algorithm: go throught the list of original faces and
    // distribute them to the cells, skipping the ones marked as
    // removed.  Then add the new and modified faces to the cells.
    // Renumber the faces using the point renumbering map.
    // Gather the internal faces by looping through the cell list and
    // collecting faces

    const labelList& allOwn = allOwner();
    const labelList& allNei = allNeighbour();

    List<DynamicList<face, primitiveMesh::facesPerCell_, 1> >
        cf(cells_.size() + ref.addedCells().size());

    List<DynamicList<label, primitiveMesh::facesPerCell_, 1> >
        cfLabels(cells_.size() + ref.addedCells().size());

    // Insert untouched internal faces
    for(label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        if (!removedFaces.found(faceI))
        {
            cf[allOwn[faceI]].append(faces_[faceI]);
            cfLabels[allOwn[faceI]].append(faceI);

            cf[allNei[faceI]].append(faces_[faceI]);
            cfLabels[allNei[faceI]].append(faceI);
        }
    }

    if (morphDebug)
    {
        Info << "Inserted untouched faces into cells" << endl;
    }

    // Add the modified internal faces
    const DynamicList<polyModifyFace>& mf = ref.modifiedFaces();

    forAll (mf, mfI)
    {
        if (!mf[mfI].isInPatch() && !mf[mfI].onlyInZone())
        {
            if (morphDebug)
            {
                // Check if the internal face is defined properly
                if
                (
                    min(mf[mfI].owner(), mf[mfI].neighbour()) < 0
                 || max(mf[mfI].owner(), mf[mfI].neighbour()) >= cf.size()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Invalid modified face " << mf[mfI].faceID()
                        << ".  Declared as internal but owner or neighbour "
                        << "are invalid." << nl
                        << "Owner: " << mf[mfI].owner()
                        << " Neighbour: " << mf[mfI].neighbour()
                        << " Max number of cells: " << cf.size()
                        << abort(FatalError);
                }
            }

            // Grab face and face label
            cf[mf[mfI].owner()].append(mf[mfI].newFace());
            cfLabels[mf[mfI].owner()].append(mf[mfI].faceID());

            cf[mf[mfI].neighbour()].append(mf[mfI].newFace());
            cfLabels[mf[mfI].neighbour()].append(mf[mfI].faceID());
        }
    }

    if (morphDebug)
    {
        Info << "Inserted modified faces into cells" << endl;
    }

    // Add the new faces
    const DynamicList<polyAddFace>& af = ref.addedFaces();

    forAll (af, afI)
    {
        if (!af[afI].isInPatch() && !af[afI].onlyInZone())
        {
            if (morphDebug)
            {
                // Check if the internal face is defined properly
                if
                (
                    min(af[afI].owner(), af[afI].neighbour()) < 0
                 || max(af[afI].owner(), af[afI].neighbour()) >= cf.size()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Invalid added face " << faces_.size() + afI
                        << ".  Declared as internal but owner or neighbour "
                        << "are invalid." << nl
                        << "Owner: " << af[afI].owner()
                        << " Neighbour: " << af[afI].neighbour()
                        << " Max number of cells: " << cf.size()
                        << abort(FatalError);
                }
            }

            // Grab face and face label
            cf[af[afI].owner()].append(af[afI].newFace());
            cfLabels[af[afI].owner()].append(faces_.size() + afI);

            cf[af[afI].neighbour()].append(af[afI].newFace());
            cfLabels[af[afI].neighbour()].append(faces_.size() + afI);
        }
    }

    if (morphDebug)
    {
        Pout<< "Inserted added faces into cells" << endl;
    }

    // All internal faces now exist on the list.  Create the new face
    // list using upper triangular search

    // Create partial point-cell addressing
    List<DynamicList<label, primitiveMesh::facesPerPoint_> > PointCells
    (
        points_.size() + ref.addedPoints().size()
    );

    boolListList usedCellFaces(cf.size());

    forAll (cf, cellI)
    {
        if (!removedCells.found(cellI))
        {
            const DynamicList<face, primitiveMesh::facesPerCell_, 1>& curFaces =
                cf[cellI];

            // Resize and reset the cell-face list at the same time
            usedCellFaces[cellI].setSize(curFaces.size());
            usedCellFaces[cellI] = false;

            // Add the cell as a neighbour to all of the points
            // of all of its faces
            forAll (curFaces, faceI)
            {
                const labelList& curFacePoints = curFaces[faceI];

                forAll (curFacePoints, pointI)
                {
                    bool found = false;

                    DynamicList<label, primitiveMesh::facesPerPoint_>& 
                        curPointCells = PointCells[curFacePoints[pointI]];

                    forAll (curPointCells, i)
                    {
                        if (curPointCells[i] == cellI)
                        {
                            // Cell already a neighbour of this point
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        curPointCells.append(cellI);
                    }
                }
            }
        }
        else
        {
            // Check if the removed cell has got any faces
            if (cf[cellI].size() > 0)
            {
                FatalErrorIn
                (
                    "void polyMesh::morph(const polyTopoChange& ref)"
                )   << "Cell " << cellI << " is marked as removed but still "
                    << "has faces.  Cell faces: " << cf[cellI]
                    << abort(FatalError);
            }
        }
    }

    List<DynamicList<label, primitiveMesh::facesPerCell_, 1> >
        newCellFaces(cf.size());

    faceList newFaces(faces_.size() + ref.faceBalance());

    // renumberFaces contains a new face label for every old and added face
    labelList renumberFaces(faces_.size() + ref.addedFaces().size(), -1);

    // faceMap contains the old face label for every new face.  For
    // inserted faces, the label will be -1.  Use with care!
    labelList faceMap(faces_.size() + ref.faceBalance(), -1);
    label nNewFaces = 0;

    // Create ordered list of faces

    // Note:
    // Insertion cannot be done in one go as the faces need to be
    // added into the list in the increasing order of neighbour
    // cells.  Therefore, all neighbours will be detected first
    // and then added in the correct order.  
    // Watch out.  Subtly different from createPolyMesh!

    forAll (cf, cellI)
    {
        const DynamicList<face, primitiveMesh::facesPerCell_, 1>&
            curFaces = cf[cellI];
        labelList neiCells(curFaces.size(), -1);
        label nNeighbours = 0;

        // For all faces ...
        forAll(curFaces, faceI)
        {
            // Skip faces that have already been matched
            if (usedCellFaces[cellI][faceI]) continue;

            bool found = false;

            const face& curFace = curFaces[faceI];

            // Get the list of labels
            const labelList& curPoints = curFace;

            // For all points
            forAll(curPoints, pointI)
            {
                // Get the list of cells sharing this point
                const DynamicList<label, primitiveMesh::facesPerPoint_>&
                    curNeighbours = PointCells[curPoints[pointI]];

                // For all neighbours
                forAll(curNeighbours, neiI)
                {
                    label curNei = curNeighbours[neiI];

                    // Reject neighbours with the lower label
                    if (curNei > cellI)
                    {
                        // Get the list of search faces
                        const DynamicList<face, primitiveMesh::facesPerCell_, 1>&
                            searchFaces = cf[curNei];

                        forAll(searchFaces, neiFaceI)
                        {
                            if (searchFaces[neiFaceI] == curFace)
                            {
                                // Match!!
                                found = true;

                                // Record the neighbour cell and face
                                neiCells[faceI] = curNei;
                                nNeighbours++;

                                break;
                            }
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            } // End of current points
        } // End of current faces

        // Add the faces in the increasing order of neighbours
        for
        (
            label neiSearch = 0;
            neiSearch < nNeighbours;
            neiSearch++
        )
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = cf.size();

            forAll (neiCells, ncI)
            {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                // Add the face to the list of faces
                newFaces[nNewFaces] = curFaces[nextNei];

                // Set cell-face and cell-neighbour-face to current face label
                newCellFaces[cellI].append(nNewFaces);
                newCellFaces[minNei].append(nNewFaces);

                // Grab the renumbering index
                renumberFaces[cfLabels[cellI][nextNei]] = nNewFaces;

                if (cfLabels[cellI][nextNei] < faces_.size())
                {
                    faceMap[nNewFaces] =
                        cfLabels[cellI][nextNei];
                }

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                // Increment number of faces counter
                nNewFaces++;
            }
            else
            {
                FatalErrorIn
                (
                    "void polyMesh::morph(const polyTopoChange& ref)"
                )   << "Error in internal face insertion"
                    << abort(FatalError);
            }
        }
    }

    label debugFaceCounter = 0;

    if (morphDebug)
    {
        Pout<< "Added internal faces.  Face count = "
            << nNewFaces << endl;

        debugFaceCounter = nNewFaces;
    }

    // Redistribute modified and newly added patch faces per patch
    List<DynamicList<label, 10> > patchModifiedFaces(boundaryMesh().size());
    List<DynamicList<label, 10> > patchAddedFaces(boundaryMesh().size());

    forAll (mf, mfI)
    {
        if (mf[mfI].isInPatch())
        {
            patchModifiedFaces[mf[mfI].patchID()].append(mfI);
        }
    }

    forAll (af, afI)
    {
        if (af[afI].isInPatch())
        {
            patchAddedFaces[af[afI].patchID()].append(afI);
        }
    }

    // For every patch add original patch faces which have not been removed
    // Remember the new patch starts and sizes
    labelList patchSizes(boundaryMesh().size(), 0);
    labelList patchStarts(boundaryMesh().size(), -1);

    forAll (boundaryMesh(), patchI)
    {
        // Add original patch faces
        const label curSize = boundaryMesh()[patchI].size();
        const label curStart = boundaryMesh()[patchI].start();

        const unallocFaceList& curFaces = boundaryMesh()[patchI];

        // Grab patch start
        patchStarts[patchI] = nNewFaces;

        for (label faceI = curStart; faceI < curStart + curSize; faceI++)
        {
            if (!removedFaces.found(faceI))
            {
                newFaces[nNewFaces] = curFaces[faceI - curStart];
                renumberFaces[faceI] = nNewFaces;

                faceMap[nNewFaces] = faceI;

                // Add the face to the owner cell
                if (allOwn[faceI] >= 0)
                {
                    newCellFaces[allOwn[faceI]].append(nNewFaces);
                }

                nNewFaces++;
            }
        }

        if (morphDebug)
        {
            Pout<< "Patch " << patchI
                << ": added faces: untouched = "
                << nNewFaces - debugFaceCounter;

            debugFaceCounter = nNewFaces;
        }

        // Add new faces belonging to this patch
        const DynamicList<label, 10>& modPatchFaces =
            patchModifiedFaces[patchI];

        forAll (modPatchFaces, faceI)
        {
            newFaces[nNewFaces] = mf[modPatchFaces[faceI]].newFace();
            renumberFaces[mf[modPatchFaces[faceI]].faceID()] = nNewFaces;

            faceMap[nNewFaces]= mf[modPatchFaces[faceI]].faceID();

            label newOwner = mf[modPatchFaces[faceI]].owner();

            // Add the face to the owner cell
            if (newOwner >= 0)
            {
                newCellFaces[newOwner].append(nNewFaces);
            }

            nNewFaces++;
        }

        if (morphDebug)
        {
            Pout<< " modified = " << nNewFaces - debugFaceCounter;

            debugFaceCounter = nNewFaces;
        }

        // Add new faces belonging to this patch
        const DynamicList<label, 10>& newPatchFaces = patchAddedFaces[patchI];

        forAll (newPatchFaces, faceI)
        {
            newFaces[nNewFaces] = af[newPatchFaces[faceI]].newFace();
            renumberFaces[faces_.size() + newPatchFaces[faceI]] = nNewFaces;

            label newOwner = af[newPatchFaces[faceI]].owner();

            // Add the face to the owner cell
            if (newOwner >= 0)
            {
                newCellFaces[newOwner].append(nNewFaces);
            }

            nNewFaces++;
        }

        if (morphDebug)
        {
            Pout<< " added = " << nNewFaces - debugFaceCounter
                << ".  Face count = " << nNewFaces << endl;

            debugFaceCounter = nNewFaces;
        }

        // Grab the new patch size
        patchSizes[patchI] = nNewFaces - patchStarts[patchI];
    }

    // Add freely-standing faces to the back of the list
    forAll (mf, mfI)
    {
        if (mf[mfI].onlyInZone())
        {
            newFaces[nNewFaces] = mf[mfI].newFace();
            renumberFaces[mf[mfI].faceID()] = nNewFaces;
            faceMap[nNewFaces]= mf[mfI].faceID();
            nNewFaces++;
        }
    }

    if (morphDebug)
    {
        Pout<< "Added zone-only faces: modified = "
            << nNewFaces - debugFaceCounter;

            debugFaceCounter = nNewFaces;
    }

    forAll (af, afI)
    {
        if (af[afI].onlyInZone())
        {
            newFaces[nNewFaces] = af[afI].newFace();
            renumberFaces[faces_.size() + afI] = nNewFaces;
            nNewFaces++;
        }
    }

    if (morphDebug)
    {
        Pout<< " added = " << nNewFaces - debugFaceCounter
            << ".  Final face count = "
            << nNewFaces << nl << endl;

            debugFaceCounter = nNewFaces;
    }

    if (debug)
    {
        if (nNewFaces != faces_.size() + ref.faceBalance())
        {
            FatalErrorIn
            (
                "void polyMesh::morph(const polyTopoChange& ref)"
            )   << "Error in face insertion.  Number of inserted faces: "
                << nNewFaces << ".  Expected "
                << faces_.size() + ref.faceBalance()
                << " faces."
                << abort(FatalError);
        }
    }

    // Face list and cell faces completed.
    // Renumber the faces using the face renumber list
    forAll (newFaces, faceI)
    {
        face oldFace = newFaces[faceI];
        face& renumberedFace = newFaces[faceI];

        forAll (renumberedFace, pointI)
        {
            renumberedFace[pointI] = renumberPoints[oldFace[pointI]];
        }

        if (morphDebug)
        {
            // Check if the face has been mapped correctly
            if
            (
                renumberedFace.size() == 0
             || min(renumberedFace) < 0
             || max(renumberedFace) >= newPointsZeroVol.size()
            )
            {
                FatalErrorIn
                (
                    "void polyMesh::morph(const polyTopoChange& ref)"
                )   << "Face " << faceI << " in the new mesh is not "
                    << "mapped correctly." << nl
                    << "It uses a removed or a non-existing vertex or "
                    << "has been skipped ." << nl
                    << "Face before mapping: " << oldFace << nl
                    << "Face after mapping: " << renumberedFace << nl
                    << "Max new vertex index: "
                    << newPointsZeroVol.size() - 1 << "." << nl
                    << "Are there extra faces in the face list that do not "
                    << "belong to a face zone?  This is not allowed."
                    << abort(FatalError);
            }
        }
    }

    // Build the face-from maps

    List<objectMap> faceFromPoint(af.size());
    label nFaceFromPoint = 0;
    List<objectMap> faceFromEdge(af.size());
    label nFaceFromEdge = 0;

    forAll (af, afI)
    {
        if (af[afI].isPointMaster())
        {
            if (morphDebug)
            {
                // Check that the master point index is in range
                if
                (
                    af[afI].masterPointID() < 0
                 || af[afI].masterPointID() >= nPoints()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Master point for face " << faces_.size() + afI
                        << " is out of range: " << af[afI].masterPointID()
                        << ".\n  Number of valid master points: "
                        << nPoints()
                        << abort(FatalError);
                }
            }

            if (af[afI].isInPatch())
            {
                // Grab faces neighbouring the point which are in the
                // same patch as the newly added face.
                const labelList& pf =
                   pointFaces()[af[afI].masterPointID()];

                labelList facesAroundPoint(pf.size());
                label nfap = 0;

                forAll (pf, pfI)
                {
                    label wp = boundaryMesh().whichPatch(pf[pfI]);
                    if (wp == af[afI].patchID())
                    {
                        facesAroundPoint[nfap] = pf[pfI];
                        nfap++;
                    }
                }

                if (morphDebug)
                {
                    if (nfap == 0)
                    {
                        FatalErrorIn
                        (
                            "void polyMesh::morph(const polyTopoChange&)"
                        )   << "No patch face neighbours found for added "
                            << "patch face " << afI
                            << ".\nThere are no faces from patch "
                            << af[afI].patchID()
                            << " around the master point "
                            << af[afI].masterPointID()
                            << ".\n Bad choice of master point: "
                            << "error in mesh mapping."
                            << abort(FatalError);
                    }
                }

                facesAroundPoint.setSize(nfap);

                faceFromPoint[nFaceFromPoint] =
                    objectMap
                    (
                        renumberFaces[faces_.size() + afI],
                        facesAroundPoint
                    );

                nFaceFromPoint++;
            }
            else
            {
                // Grab internal faces around the point
                const labelList& pf =
                    pointFaces()[af[afI].masterPointID()];

                labelList facesAroundPoint(pf.size());
                label nfap = 0;

                forAll (pf, pfI)
                {
                    if (isInternalFace(pf[pfI]))
                    {
                        facesAroundPoint[nfap] = pf[pfI];
                        nfap++;
                    }
                }

                if (morphDebug)
                {
                    if (nfap == 0 && nInternalFaces() > 0)
                    {
                        FatalErrorIn
                        (
                            "void polyMesh::morph(const polyTopoChange&)"
                        )   << "No face neighbours found for added "
                            << "internal face " << afI
                            << ".\nThere are no internal faces "
                            << "around the master point "
                            << af[afI].masterPointID()
                            << ".\n Bad choice of master point: "
                            << "error in mesh mapping."
                            << abort(FatalError);
                    }
                }

                facesAroundPoint.setSize(nfap);

                faceFromPoint[nFaceFromPoint] =
                    objectMap
                    (
                        renumberFaces[faces_.size() + afI],
                        facesAroundPoint
                    );

                nFaceFromPoint++;
            }
        }
        else if (af[afI].isEdgeMaster())
        {
            if (morphDebug)
            {
                // Check that the master edge index is in range
                if
                (
                    af[afI].masterEdgeID() < 0
                 || af[afI].masterEdgeID() >= nEdges()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Master edge for face " << faces_.size() + afI
                        << " is out of range: " << af[afI].masterEdgeID()
                        << ".\n  Number of valid master edges: "
                        << nEdges()
                        << abort(FatalError);
                }
            }

            if (af[afI].isInPatch())
            {
                // Grab faces neighbouring the point which are in the
                // same patch as the newly added face Note: the
                // addressing is now into the patch instead of the
                // global face list
                const labelList& pe =
                   edgeFaces()[af[afI].masterEdgeID()];

                labelList facesAroundEdge(pe.size());
                label nfae = 0;

                forAll (pe, peI)
                {
                    label wp = boundaryMesh().whichPatch(pe[peI]);
                    if (wp == af[afI].patchID())
                    {
                        facesAroundEdge[nfae] = pe[peI];
                        nfae++;
                    }
                }

                if (morphDebug)
                {
                    if (nfae == 0)
                    {
                        FatalErrorIn
                        (
                            "void polyMesh::morph(const polyTopoChange&)"
                        )   << "No patch face neighbours found for added "
                            << "patch face " << afI
                            << ".  Error in mesh mapping."
                            << abort(FatalError);
                    }
                }

                facesAroundEdge.setSize(nfae);

                faceFromEdge[nFaceFromEdge] =
                    objectMap
                    (
                        renumberFaces[faces_.size() + afI],
                        facesAroundEdge
                    );

                nFaceFromEdge++;
            }
            else
            {
                // Grab internal faces around the edge
                const labelList& pe =
                    edgeFaces()[af[afI].masterEdgeID()];

                labelList facesAroundEdge(pe.size());
                label nfae = 0;

                forAll (pe, peI)
                {
                    if (isInternalFace(pe[peI]))
                    {
                        facesAroundEdge[nfae] = pe[peI];
                        nfae++;
                    }
                }

                if (morphDebug)
                {
                    if (nfae == 0)
                    {
                        FatalErrorIn
                        (
                            "void polyMesh::morph(const polyTopoChange&)"
                        )   << "No patch face neighbours found for added "
                            << "internal face " << afI
                            << ".  Error in mesh mapping."
                            << abort(FatalError);
                    }
                }

                facesAroundEdge.setSize(nfae);

                faceFromEdge[nFaceFromEdge] =
                    objectMap
                    (
                        renumberFaces[faces_.size() + afI],
                        facesAroundEdge
                    );

                nFaceFromEdge++;
            }
        }
        else if (af[afI].isFaceMaster()) // Face mastered by another face
        {
            faceMap[renumberFaces[faces_.size() + afI]] =
                af[afI].masterFaceID();
        }
    }

    // Reset the size of face mapping lists
    faceFromPoint.setSize(nFaceFromPoint);
    faceFromEdge.setSize(nFaceFromEdge);

    // Check face maps
    if (morphDebug)
    {
        boolList mappedFaces(faceMap.size(), false);

        // Fill in faces mapped from the face map
        forAll (faceMap, faceI)
        {
            if (faceMap[faceI] >= 0)
            {
                mappedFaces[faceI] = true;
            }
        }

        // Fill in point and edge maps
        forAll (faceFromPoint, faceI)
        {
            mappedFaces[faceFromPoint[faceI].index()] = true;
        }

        forAll (faceFromEdge, faceI)
        {
            mappedFaces[faceFromEdge[faceI].index()] = true;
        }

        // Check if all the faces are mapped
        label nUnmappedFaces = 0;

        forAll (mappedFaces, faceI)
        {
            if (!mappedFaces[faceI])
            {
                nUnmappedFaces++;
            }
        }

        if (nUnmappedFaces > 0)
        {
            Pout<< "void polyMesh::morph(const polyTopoChange& ref) : "
                << "unmapped data for " << nUnmappedFaces << " faces." << endl;
        }
    }

    labelHashSet flipFaceFlux(mf.size() + af.size());

    // Build the flip flux map
    forAll (mf, mfI)
    {
        if (mf[mfI].flipFaceFlux())
        {
            flipFaceFlux.insert(renumberFaces[mf[mfI].faceID()]);
        }
    }

    forAll (af, afI)
    {
        if (af[afI].isFaceMaster() && af[afI].flipFaceFlux())
        {
            flipFaceFlux.insert(renumberFaces[faces_.size() + afI]);
        }
    }

    // Renumber the cells

    cellList newCells(cells_.size() + ref.cellBalance());

    // renumberCells holds the new cell label for all old and added cells
    labelList renumberCells(cells_.size() + ref.addedCells().size(), -1);

    // cellMap holds the old cell label for every preserved cell
    labelList cellMap(cells_.size() + ref.cellBalance(), -1);
    
    label nNewCells = 0;

    forAll (newCellFaces, cellI)
    {
        if (!removedCells.found(cellI))
        {
            if (newCellFaces[cellI].size() < 4)
            {
                FatalErrorIn
                (
                    "void polyMesh::morph(const polyTopoChange& ref)"
                )   << "Cell " << cellI << " has got three or less faces "
                    << "and has not been removed.  "
                    << "This is not a valid cell." << endl
                    << "Cell faces: " << newCellFaces[cellI]
                    << abort(FatalError);
            }

            // Make a cell
            newCells[nNewCells].transfer(newCellFaces[cellI].shrink());

            renumberCells[cellI] = nNewCells;

            // Add the cell into cell map if it is preserved
            if (cellI < nOldCells)
            {
                cellMap[nNewCells] = cellI;
            }

            nNewCells++;
        }
    }

    if (morphDebug)
    {
        Pout<< "Added all cells.  Final cell count = "
            << nNewCells << nl << endl;
    }

    label nPreservedCells = cells_.size() - removedCells.size();

    // Build the cell-from maps

    const DynamicList<polyAddCell>& ac = ref.addedCells();
    const DynamicList<polyModifyCell>& mc = ref.modifiedCells();

    List<objectMap> cellFromPoint(ac.size());
    label nCellFromPoint = 0;

    List<objectMap> cellFromEdge(ac.size());
    label nCellFromEdge = 0;

    List<objectMap> cellFromFace(ac.size());
    label nCellFromFace = 0;

    forAll (ac, acI)
    {
        if (ac[acI].isPointMaster())
        {
            if (morphDebug)
            {
                // Check that the master point index is in range
                if
                (
                    ac[acI].masterPointID() < 0
                 || ac[acI].masterPointID() >= nPoints()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Master point for cell " << nPreservedCells + acI
                        << " is out of range: " << ac[acI].masterPointID()
                        << ".\n  Number of valid master points: "
                        << nPoints()
                        << abort(FatalError);
                }
            }

            cellFromPoint[nCellFromPoint] =
                objectMap
                (
                    nPreservedCells + acI,
                    pointCells()[ac[acI].masterPointID()]
                );

            nCellFromPoint++;
        }
        else if (ac[acI].isEdgeMaster())
        {
            if (morphDebug)
            {
                // Check that the master edge index is in range
                if
                (
                    ac[acI].masterEdgeID() < 0
                 || ac[acI].masterEdgeID() >= nEdges()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Master edge for cell " << nPreservedCells + acI
                        << " is out of range: " << ac[acI].masterEdgeID()
                        << ".\n  Number of valid master edges: "
                        << nEdges()
                        << abort(FatalError);
                }
            }

            cellFromEdge[nCellFromEdge] =
                objectMap
                (
                    nPreservedCells + acI,
                    edgeCells()[ac[acI].masterEdgeID()]
                );

            nCellFromEdge++;
        }
        else if (ac[acI].isFaceMaster()) // Cell mastered by a face
        {
            if (morphDebug)
            {
                // Check that the master face index is in range
                if
                (
                    ac[acI].masterFaceID() < 0
                 || ac[acI].masterFaceID() >= nFaces()
                )
                {
                    FatalErrorIn
                    (
                        "void polyMesh::morph(const polyTopoChange& ref)"
                    )   << "Master face for cell " << nPreservedCells + acI
                        << " is out of range: " << ac[acI].masterFaceID()
                        << ".\n  Number of valid master faces: "
                        << nFaces()
                        << abort(FatalError);
                }
            }

            labelList cellsAroundFace(2, -1);

            cellsAroundFace[0] = allOwner()[ac[acI].masterFaceID()];

            if (allNeighbour()[ac[acI].masterFaceID()] >= 0)
            {
                cellsAroundFace[1] = allNeighbour()[ac[acI].masterFaceID()];
            }
            else
            {
                cellsAroundFace.setSize(1);
            }

            cellFromFace[nCellFromFace] =
                objectMap
                (
                    nPreservedCells + acI,
                    cellsAroundFace
                );

            nCellFromFace++;
        }
        else if (ac[acI].isCellMaster()) // Cell mastered by another cell
        {
            cellMap[renumberCells[cells_.size() + acI]] =
                ac[acI].masterCellID();
        }
    }

    // Reset the size of cell mapping lists
    cellFromPoint.setSize(nCellFromPoint);
    cellFromEdge.setSize(nCellFromEdge);
    cellFromFace.setSize(nCellFromFace);

    // Check cell maps
    if (morphDebug)
    {
        boolList mappedCells(newCells.size(), false);

        // Fill in cells mapped from the cell map
        forAll (cellMap, cellI)
        {
            if (cellMap[cellI] >= 0)
            {
                mappedCells[cellI] = true;
            }
        }

        // Fill in point and edge maps
        forAll (cellFromPoint, cellI)
        {
            mappedCells[cellFromPoint[cellI].index()] = true;
        }

        forAll (cellFromEdge, cellI)
        {
            mappedCells[cellFromEdge[cellI].index()] = true;
        }

        forAll (cellFromFace, cellI)
        {
            mappedCells[cellFromFace[cellI].index()] = true;
        }

        // Check if all the cells are mapped
        label nUnmappedCells = 0;

        forAll (mappedCells, cellI)
        {
            if (!mappedCells[cellI])
            {
                nUnmappedCells++;
            }
        }

        if (nUnmappedCells > 0)
        {
            Pout<< "void polyMesh::morph(const polyTopoChange& ref) : "
                << "unmapped data for " << nUnmappedCells << " cells." << endl;
        }
    }


    // Rebuild the mesh
    // ~~~~~~~~~~~~~~~~

    // Grab patch mesh point maps

    List<Map<label> > oldPatchMeshPointMaps(boundary_.size());
    labelList oldPatchNMeshPoints(boundary_.size());

    forAll (boundary_, patchI)
    {
        // Copy old face zone mesh point maps
        oldPatchMeshPointMaps[patchI] = boundary_[patchI].meshPointMap();
        oldPatchNMeshPoints[patchI] = boundary_[patchI].meshPoints().size();
    }

    // Re-do the point zones

    // Make a map of points to be removed from zones
    labelHashSet removePointFromZone(2*ref.modifiedPoints().size());

    forAll (mp, mpI)
    {
        if (mp[mpI].removeFromZone())
        {
            removePointFromZone.insert(mp[mpI].pointID());
        }
    }

    labelListList newPointZoneAddr(pointZones_.size());
    labelList nPointsInZone(pointZones_.size(), 0);

    forAll (pointZones_, pzI)
    {
        // Get the list of old points
        const labelList& oldAddr = pointZones_[pzI].addressing();

        // Create new addressing, over-estimating the size
        labelList& newAddr = newPointZoneAddr[pzI];
        label& curNPoints = nPointsInZone[pzI];

        newAddr.setSize
        (
            oldAddr.size()
          + ref.modifiedPoints().size()
          + ref.addedPoints().size()
        );

        // Add the original points that have not been removed or re-zoned
        forAll (oldAddr, pointI)
        {
            if
            (
                !removedPoints.found(oldAddr[pointI])
             && !removePointFromZone.found(oldAddr[pointI])
            )
            {
                // The point is still alive. Add its renumbered label
                newAddr[curNPoints] = renumberPoints[oldAddr[pointI]];
                curNPoints++;
            }
        }
    }

    labelList debugPointsInZone(pointZones_.size(), 0);

    if (morphDebug)
    {
        Pout<< "Added zone points:  untouched = " << nPointsInZone;

        debugPointsInZone = nPointsInZone;
    }

    // Distribute modified zone points
    forAll (mp, mpI)
    {
        if (mp[mpI].isInZone())
        {
            const label zoneID = mp[mpI].zoneID();

            newPointZoneAddr[zoneID][nPointsInZone[zoneID]] =
                renumberPoints[mp[mpI].pointID()];

            nPointsInZone[zoneID]++;
        }
    }

    if (morphDebug)
    {
        Pout<< " modified = " << nPointsInZone - debugPointsInZone;
    }

    // Distribute added zone points
    forAll (ap, apI)
    {
        if (ap[apI].isInZone())
        {
            const label zoneID = ap[apI].zoneID();

            newPointZoneAddr[zoneID][nPointsInZone[zoneID]] =
                renumberPoints[points_.size() + apI];

            nPointsInZone[zoneID]++;
        }
    }

    if (morphDebug)
    {
        Pout<< " added = " << nPointsInZone - debugPointsInZone
            << ".  Points per zone = " << nPointsInZone << endl;
    }

    // Reset the sizes of the point zone addressing
    forAll (newPointZoneAddr, pzI)
    {
        newPointZoneAddr[pzI].setSize(nPointsInZone[pzI]);
    }

    // Build the point zone renumbering
    labelListList pzRenumber(pointZones_.size());

    forAll (pointZones_, pzI)
    {
        pointZone& oldZone = pointZones_[pzI];
        const labelList& newZoneAddr = newPointZoneAddr[pzI];

        labelList& curPzRnb = pzRenumber[pzI];
        curPzRnb.setSize(newZoneAddr.size());

        forAll (newZoneAddr, pointI)
        {
            if (newZoneAddr[pointI] < pointMap.size())
            {
                curPzRnb[pointI] =
                    oldZone.whichPoint(pointMap[newZoneAddr[pointI]]);
            }
            else
            {
                curPzRnb[pointI] = -1;
            }
        }
    }


    // Re-do the face zones

    // Make a map of faces to be removed from zones
    labelHashSet removeFaceFromZone(2*ref.modifiedFaces().size());

    forAll (mf, mfI)
    {
        if (mf[mfI].removeFromZone())
        {
            removeFaceFromZone.insert(mf[mfI].faceID());
        }
    }

    labelListList newFaceZoneAddr(faceZones_.size());
    boolListList newFaceZoneFaceFlip(faceZones_.size());
    labelList nFacesInZone(faceZones_.size(), 0);

    forAll (faceZones_, fzI)
    {
        // Get the list of old faces
        const labelList& oldAddr = faceZones_[fzI].addressing();
        const boolList& oldFlip = faceZones_[fzI].flipMap();

        // Create new addressing, over-estimating the size
        labelList& newAddr = newFaceZoneAddr[fzI];
        boolList& newFlip = newFaceZoneFaceFlip[fzI];
        label& curNFaces = nFacesInZone[fzI];

        newAddr.setSize
        (
            oldAddr.size()
          + ref.modifiedFaces().size()
          + ref.addedFaces().size()
        );

        newFlip.setSize(newAddr.size());
        newFlip = false;

        // Reset the face flip to false.  None of the preserved faces will
        // be flipped

        // Add the original faces that have not been removed or re-zoned
        forAll (oldAddr, faceI)
        {
            if
            (
                !removedFaces.found(oldAddr[faceI])
             && !removeFaceFromZone.found(oldAddr[faceI])
            )
            {
                // The face is still alive. Add its renumbered label
                newAddr[curNFaces] = renumberFaces[oldAddr[faceI]];
                newFlip[curNFaces] = oldFlip[faceI];

                curNFaces++;
            }
        }
    }

    labelList debugFacesInZone(faceZones_.size(), 0);

    if (morphDebug)
    {
        Pout<< "Added zone faces: untouched = " << nFacesInZone;

        debugFacesInZone = nFacesInZone;
    }

    // Distribute modified zone faces
    forAll (mf, mfI)
    {
        if (mf[mfI].isInZone())
        {
            const label zoneID = mf[mfI].zoneID();

            // Grab the face index
            newFaceZoneAddr[zoneID][nFacesInZone[zoneID]] =
                renumberFaces[mf[mfI].faceID()];

            // Grab the face flip
            newFaceZoneFaceFlip[zoneID][nFacesInZone[zoneID]] =
                mf[mfI].zoneFlip();

            nFacesInZone[zoneID]++;
        }
    }

    if (morphDebug)
    {
        Pout<< " modified = " << nFacesInZone - debugFacesInZone;

        debugFacesInZone = nFacesInZone;
    }

    // Distribute added zone faces
    forAll (af, afI)
    {
        if (af[afI].isInZone())
        {
            const label zoneID = af[afI].zoneID();

            // Grab the face index
            newFaceZoneAddr[zoneID][nFacesInZone[zoneID]] =
                renumberFaces[faces_.size() + afI];

            newFaceZoneFaceFlip[zoneID][nFacesInZone[zoneID]] =
                af[afI].zoneFlip();

            nFacesInZone[zoneID]++;
        }
    }

    if (morphDebug)
    {
        Pout<< " added = " << nFacesInZone - debugFacesInZone
            << ".  Faces per zone = " << nFacesInZone << endl;
    }

    // Reset the sizes of the face zone addressing and face flip
    forAll (newFaceZoneAddr, fzI)
    {
        newFaceZoneAddr[fzI].setSize(nFacesInZone[fzI]);
        newFaceZoneFaceFlip[fzI].setSize(nFacesInZone[fzI]);
    }

    // Build the face zone renumbering
    labelListList fzFaceRenumber(faceZones_.size());

    List<Map<label> > oldFaceZoneMeshPointMaps(faceZones_.size());

    forAll (faceZones_, fzI)
    {
        faceZone& oldZone = faceZones_[fzI];

        // Point renumbering
        // Point renumbering gives the old point location of every new
        // point in the face zone.  If the point is new to the zone,
        // the index will be -1.
        // The problem is that the renumbering cannot be built at this stage
        // as the order of points in the new face zone is not yet known:
        // (the new zone needs to be hooked onto the updated list of faces to
        // be able to create the mapping).  Therefore, the old meshPoint maps
        // will be copied here and will be used later to re-create the
        // addressing.  

        // Copy old face zone mesh point maps
        oldFaceZoneMeshPointMaps[fzI] = faceZones_[fzI]().meshPointMap();

        // Face renumbering
        const labelList& newZoneAddr = newFaceZoneAddr[fzI];

        labelList& curFzFaceRnb = fzFaceRenumber[fzI];

        curFzFaceRnb.setSize(newZoneAddr.size());

        forAll (newZoneAddr, faceI)
        {
            if (newZoneAddr[faceI] < faceMap.size())
            {
                curFzFaceRnb[faceI] =
                    oldZone.whichFace(faceMap[newZoneAddr[faceI]]);
            }
            else
            {
                curFzFaceRnb[faceI] = -1;
            }
        }
    }


    // Re-do the cell zones

    // Make a map of cells to be removed from zones
    labelHashSet removeCellFromZone(2*ref.modifiedCells().size());

    forAll (mc, mcI)
    {
        if (mc[mcI].removeFromZone())
        {
            removeCellFromZone.insert(mc[mcI].cellID());
        }
    }

    labelListList newCellZoneAddr(cellZones_.size());
    labelList nCellsInZone(cellZones_.size(), 0);

    forAll (cellZones_, czI)
    {
        // Get the list of old cells
        const labelList& oldAddr = cellZones_[czI].addressing();

        // Create new addressing, over-estimating the size
        labelList& newAddr = newCellZoneAddr[czI];
        label& curNCells = nCellsInZone[czI];

        newAddr.setSize
        (
            oldAddr.size()
          + ref.modifiedCells().size()
          + ref.addedCells().size()
        );

        // Add the original cells that have not been removed or re-zoned
        forAll (oldAddr, cellI)
        {
            if
            (
                !removedCells.found(oldAddr[cellI])
             && !removeCellFromZone.found(oldAddr[cellI])
            )
            {
                // The cell is still alive. Add its renumbered label
                newAddr[curNCells] = renumberCells[oldAddr[cellI]];
                curNCells++;
            }
        }
    }

    labelList debugCellsInZone(cellZones_.size(), 0);

    if (morphDebug)
    {
        Pout<< "Added zone cells: untouched = " << nCellsInZone;

        debugCellsInZone = nCellsInZone;
    }

    // Distribute modified zone cells
    forAll (mc, mcI)
    {
        if (mc[mcI].isInZone())
        {
            const label zoneID = mc[mcI].zoneID();

            newCellZoneAddr[zoneID][nCellsInZone[zoneID]] =
                renumberCells[mc[mcI].cellID()];

            nCellsInZone[zoneID]++;
        }
    }

    if (morphDebug)
    {
        Pout<< " modified = " << nCellsInZone - debugCellsInZone;

        debugCellsInZone = nCellsInZone;
    }

    // Distribute added zone cells
    forAll (ac, acI)
    {
        if (ac[acI].isInZone())
        {
            const label zoneID = ac[acI].zoneID();

            newCellZoneAddr[zoneID][nCellsInZone[zoneID]] =
                renumberCells[cells_.size() + acI];

            nCellsInZone[zoneID]++;
        }
    }

    if (morphDebug)
    {
        Pout<< " added = " << nCellsInZone - debugCellsInZone
            << ".  Cells per zone = " << nCellsInZone << endl;
    }

    // Reset the sizes of the cell zone addressing
    forAll (newCellZoneAddr, czI)
    {
        newCellZoneAddr[czI].setSize(nCellsInZone[czI]);
    }

    // Build the cell zone renumbering
    labelListList czRenumber(cellZones_.size());

    forAll (cellZones_, czI)
    {
        cellZone& oldZone = cellZones_[czI];
        const labelList& newZoneAddr = newCellZoneAddr[czI];

        labelList& curCzRnb = czRenumber[czI];

        curCzRnb.setSize(newZoneAddr.size());

        forAll (newZoneAddr, cellI)
        {
            if (newZoneAddr[cellI] < cellMap.size())
            {
                curCzRnb[cellI] =
                    oldZone.whichCell(cellMap[newZoneAddr[cellI]]);
            }
            else
            {
                curCzRnb[cellI] = -1;
            }
        }
    }


    // Reset the new points, faces, cells
    points_ = newPointsZeroVol;
    faces_ = newFaces;
    cells_ = newCells;

    // Reset the boundary patches
    forAll (boundary_, patchI)
    {
        boundary_[patchI] = polyPatch
        (
            boundary_[patchI].name(),
            patchSizes[patchI],
            patchStarts[patchI],
            patchI,
            boundary_
        );
    }

    // Recalculate the owner/neighbour addressing and reset the primitiveMesh
    clearFaceCells();
    calcFaceCells();

    // Reset the zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    forAll (pointZones_, pzI)
    {
        pointZones_[pzI].resetAddressing(newPointZoneAddr[pzI]);
    }

    forAll (faceZones_, fzI)
    {
        faceZones_[fzI].resetAddressing
        (
            newFaceZoneAddr[fzI],
            newFaceZoneFaceFlip[fzI]
        );
    }

    forAll (cellZones_, czI)
    {
        cellZones_[czI].resetAddressing(newCellZoneAddr[czI]);
    }


    // Map the old motion points if present
    if (oldPointsPtr_)
    {
        // Make a copy of the original points
        pointField oldMotionPoints = *oldPointsPtr_;

        pointField& newMotionPoints = *oldPointsPtr_;

        // Resize the list to new size
        newMotionPoints.setSize(points_.size());

        // Map the list
        newMotionPoints.map(oldMotionPoints, pointMap);
    }

    // Change the file options for the moving mesh

    points_.writeOpt() = IOobject::AUTO_WRITE;
    points_.instance() = time().timeName();

    faces_.writeOpt() = IOobject::AUTO_WRITE;
    faces_.instance() = time().timeName();

    cells_.writeOpt() = IOobject::AUTO_WRITE;
    cells_.instance() = time().timeName();

    boundary_.writeOpt() = IOobject::AUTO_WRITE;
    boundary_.instance() = time().timeName();

    pointZones_.writeOpt() = IOobject::AUTO_WRITE;
    pointZones_.instance() = time().timeName();

    faceZones_.writeOpt() = IOobject::AUTO_WRITE;
    faceZones_.instance() = time().timeName();

    cellZones_.writeOpt() = IOobject::AUTO_WRITE;
    cellZones_.instance() = time().timeName();


    // Create the patch mesh point renumbering

    labelListList patchPointRenumber(boundary_.size());

    forAll (boundary_, patchI)
    {
        const labelList& newPatchMeshPoints = boundary_[patchI].meshPoints();

        const Map<label>& oldZoneMeshPointMap = oldPatchMeshPointMaps[patchI];
        const label oldSize = oldPatchNMeshPoints[patchI];

        labelList& curPatchPointRnb = patchPointRenumber[patchI];

        curPatchPointRnb.setSize(newPatchMeshPoints.size());

        forAll (newPatchMeshPoints, pointI)
        {
            if (newPatchMeshPoints[pointI] < oldSize)
            {
                Map<label>::const_iterator ozmpmIter =
                    oldZoneMeshPointMap.find
                    (
                        pointMap[newPatchMeshPoints[pointI]]
                    );

                if (ozmpmIter != oldZoneMeshPointMap.end())
                {
                    curPatchPointRnb[pointI] = ozmpmIter();
                }
                else
                {
                    curPatchPointRnb[pointI] = -1;
                }
            }
            else
            {
                curPatchPointRnb[pointI] = -1;
            }
        }
    }

    // Create the face zone mesh point renumbering

    labelListList fzPointRenumber(faceZones_.size());

    forAll (faceZones_, fzI)
    {
        const labelList& newZoneMeshPoints = faceZones_[fzI]().meshPoints();

        const Map<label>& oldZoneMeshPointMap = oldFaceZoneMeshPointMaps[fzI];

        labelList& curFzPointRnb = fzPointRenumber[fzI];

        curFzPointRnb.setSize(newZoneMeshPoints.size());

        forAll (newZoneMeshPoints, pointI)
        {
            if (newZoneMeshPoints[pointI] < pointMap.size())
            {
                Map<label>::const_iterator ozmpmIter =
                    oldZoneMeshPointMap.find
                    (
                        pointMap[newZoneMeshPoints[pointI]]
                    );

                if (ozmpmIter != oldZoneMeshPointMap.end())
                {
                    curFzPointRnb[pointI] = ozmpmIter();
                }
                else
                {
                    curFzPointRnb[pointI] = -1;
                }
            }
            else
            {
                curFzPointRnb[pointI] = -1;
            }
        }
    }


    if (debug || morphDebug)
    {
        Pout<< "Foam::polyMesh::morph" << nl
            << '(' << nl
            << "    const polyTopoChange& ref" << nl
            << ") : completed topological change." << nl << endl;
    }

    morphMap_ = new mapPolyMesh
    (
        *this,
        nOldPoints,
        nOldFaces,
        nOldCells,
        pointMap,
        faceMap,
        faceFromPoint,
        faceFromEdge,
        cellMap,
        cellFromPoint,
        cellFromEdge,
        cellFromFace,
        renumberPoints,
        renumberFaces,
        renumberCells,
        flipFaceFlux,
        patchPointRenumber,
        pzRenumber,
        fzPointRenumber,
        fzFaceRenumber,
        czRenumber,
        newPointsMotion,
        oldPatchStarts,
        oldPatchNMeshPoints
    );


    if (hasCoupled)
    {
        (void)reorderCoupledPatches(ref, *morphMap_);
    }
}


// ************************************************************************* //
