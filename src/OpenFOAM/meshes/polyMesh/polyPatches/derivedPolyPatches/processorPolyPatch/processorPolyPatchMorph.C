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
    Implements patch face ordering on processor patches.
    See cyclicPolyPatchMorph as well.

    Does geometric matching (assumes face centres are coincident i.e. no cyclics
    allowed) or topological
    matching (needs at least one unmodified face per connected patch area)

    master: lowest numbered processor


    sendOrder:
        - master:
            geometric:get face centres and f[0] coords. Send to neighbour.
            topological:get seed faces and walk order. Send to neighbour.
        - slave: nothing

    order:
        - master: nothing. Ordering taken as ok already.
        - slave:
            geometric: receive coords from  master. Do geometric matching
            on coords.
            topological: receive walk order from  master. Do walk itself
            and combine the two walks.

    The difference with cyclics is that there the half0 also reorders itself
    (at least for topological matching)
    We cannot do this on processor patches since this would require the
    seed faces to be synchronized on both patches (since both of them would have
    to be unmodified). Synchronization would require an extra communication
    phase (since already in sendOrder the master would need the synchronized
    seed faces). So instead we choose to consider the ordering of the master
    patch to be ok and only do something to the slave.

\*---------------------------------------------------------------------------*/

#include "processorPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "IPstream.H"
#include "OPstream.H"
#include "walkPatch.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void processorPolyPatch::sendGeometricOrder
(
    const polyTopoChange&,
    const mapPolyMesh& map
) const
{
    const pointField& preMotionPoints = map.preMotionPoints();

    pointField ctrs(calcFaceCentres(*this, preMotionPoints));

    pointField anchorPoints(getAnchorPoints(*this, preMotionPoints));

    // Now send all info over to the neighbour
    {
        if (polyMesh::morphDebug)
        {
            Sout<< "processorPolyPatch::sendOrder :"
                << " sending face centres:" << ctrs.size()
                << " sending anchors:" << anchorPoints.size()
                << endl;
        }
        OPstream toNeighbour(neighbProcNo());
        toNeighbour << ctrs << anchorPoints;
    }
}


bool processorPolyPatch::geometricOrder
(
    const polyTopoChange&,
    const mapPolyMesh& map,
    labelList& faceMap,
    labelList& rotation
) const
{
    const pointField& preMotionPoints = map.preMotionPoints();

    vectorField masterCtrs;
    vectorField masterAnchorPoints;

    // Receive data from neighbour
    {
        IPstream fromNeighbour(neighbProcNo());
        fromNeighbour >> masterCtrs >> masterAnchorPoints;
    }

    // Calculate my face centres
    pointField ctrs(calcFaceCentres(*this, preMotionPoints));

    // Calculate typical distance from face centre
    scalarField tols(calcFaceTol(*this, preMotionPoints, ctrs));

    // Geometric match of face centre vectors
    bool matchedAll = matchPoints(ctrs, masterCtrs, tols, false, faceMap);

    if (polyMesh::morphDebug)
    {
        fileName ccName(name() + "_faceCentres.obj");

        Sout<< "processorPolyPatch::order : "
            << "Dumping newly found match as lines between"
            << " corresponding face centres to OBJ file " << ccName
            << endl;

        OFstream ccStr(ccName);

        label vertI = 0;

        forAll(ctrs, faceI)
        {
            label masterFaceI = faceMap[faceI];

            if (masterFaceI != -1)
            {
                const point& c0 = masterCtrs[masterFaceI];
                const point& c1 = ctrs[faceI];
                writeOBJ(ccStr, c0, c1, vertI);
            }
        }
    }

    if (!matchedAll)
    {
        SeriousError
            << "processorPolyPatch::order : patch:" << name() << " : "
            << "Cannot match vectors to faces on both sides of patch" << endl
            << "masterCtrs[0]:" << masterCtrs[0] << endl
            << "ctrs[0]:" << ctrs[0] << endl
            << "Please use topological ordering or adapt your topology changes"
            << endl
            << "Continuing with incorrect face ordering from now on!"
            << endl;

        return false;
    }

    // Set rotation.
    forAll(faceMap, oldFaceI)
    {
        // The face f will be at newFaceI (after morphing) and we want its
        // anchorPoint (= f[0]) to align with the anchorpoint for the
        // corresponding face on the other side.

        label newFaceI = faceMap[oldFaceI];

        const point& wantedAnchor = masterAnchorPoints[newFaceI];

        rotation[newFaceI] =
            getRotation
            (
                preMotionPoints,
                operator[](oldFaceI),
                wantedAnchor,
                tols[oldFaceI]
            );

        if (rotation[newFaceI] == -1)
        {
            SeriousError
                << "processorPolyPatch::order : patch:" << name()
                << " : "
                << "Cannot find point on face " << operator[](oldFaceI)
                << " with vertices:"
                << IndirectList<point>(preMotionPoints, operator[](oldFaceI))
                << " that matches point " << wantedAnchor
                << " when matching the halves of processor patch " << name()
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }
    }

    forAll(faceMap, faceI)
    {
        if (faceMap[faceI] != faceI || rotation[faceI] != 0)
        {
            return true;
        }
    }
    return false;
}


void processorPolyPatch::sendTopologicalOrder
(
    const polyTopoChange& ref,
    const mapPolyMesh& map
) const
{
    const labelList& oldToNew = map.reverseFaceMap();
    const labelList& oldPatchStarts = map.oldPatchStarts();

    const bool isMaster = Pstream::myProcNo() < neighbProcNo();

    // Get modified/unmodified status for every face on old mesh.
    boolList unmodifiedOldFace(getUnmodifiedFaces(ref, map));

    DynamicList<label> unmodOldPatchFaces(size());

    label oldPatchI = -1;

    forAll(unmodifiedOldFace, oldFaceI)
    {
        if (unmodifiedOldFace[oldFaceI] && inPatch(oldToNew, oldFaceI))
        {
            if (oldPatchI == -1)
            {
                // Unmodified face in the old patch (since was in
                // processorpatch)
                oldPatchI = whichPatch(oldPatchStarts, oldFaceI);

                if (oldPatchI == -1)
                {
                    FatalErrorIn
                    (
                        "processorPolyPatch::sendOrder"
                        "(const polyTopoChange&, const mapPolyMesh&)"
                    )   << "Unmodified face " << oldToNew[oldFaceI]
                        << " on processor patch " << name()
                        << " was not in a patch in old mesh"
                        << abort(FatalError);
                }                  
            }

            // Store patch in old patch relative ordering.
            unmodOldPatchFaces.append
            (
                oldFaceI
              - oldPatchStarts[oldPatchI]
            );
        }
    }
    unmodOldPatchFaces.shrink();


    // All visited faces
    boolList visited(size(), false);

    // All patch is one region
    const labelList ppZones(size(), 0);

    // All faces visited and index in face of visit.
    labelListList allVisitOrder(unmodOldPatchFaces.size());
    labelListList allIndex(unmodOldPatchFaces.size());

    forAll(unmodOldPatchFaces, seedI)
    {
        // Walk on patch starting from unmodified faces to get relative
        // connectivity

        // old mesh face label
        label oldFaceI =
            oldPatchStarts[oldPatchI] + unmodOldPatchFaces[seedI];

        // new patch face label
        label startFaceI = oldToNew[oldFaceI] - start();

        walkPatch patchWalker
        (
            *this,
            ppZones,
            !isMaster,
            startFaceI,
            localFaces()[startFaceI][0],
            visited
        );

        // Collect visiting order for this seed.
        allVisitOrder[seedI] = patchWalker.visitOrder();
        allIndex[seedI] = patchWalker.indexInFace();


        if (polyMesh::morphDebug)
        {
            if (patchWalker.visitOrder().size() > 0)
            {
                Sout<< "processorPolyPatch::sendOrder : "
                    << "On master walked from unmodified face " << startFaceI
                    << " and visited " << patchWalker.visitOrder().size()
                    << " faces." << endl;

                pointField ctrs(calcFaceCentres(*this, map.preMotionPoints()));

                fileName ccName(name() + '_' + Foam::name(seedI) + ".obj");

                Sout<< "processorPolyPatch::sendOrder : Dumping visit order for"
                    << " seed " << startFaceI << " to OBJ file " << ccName
                    << endl;

                OFstream ccStr(ccName);
                writeOBJ(ccStr, ctrs, allVisitOrder[seedI]);
            }
        }
    }

    // Now send all info over to the neighbour
    {
        if (polyMesh::morphDebug)
        {
            Sout<< "processorPolyPatch::sendOrder :"
                << " Sending seedFaces:" << unmodOldPatchFaces.size()
                << " sending visitOrder:" << allVisitOrder.size()
                << endl;
        }
        OPstream toNeighbour(neighbProcNo());
        toNeighbour << unmodOldPatchFaces << allVisitOrder << allIndex;
    }

    label falseIndex = findIndex(visited, false);

    if (falseIndex != -1)
    {
        SeriousError
            << "processorPolyPatch::sendOrder : patch:" << name()
            << " : " << "Did not visit mesh face "
            << start() + falseIndex
            << " on processor patch " << name()
            << " from any of the seed faces." << endl
            << "There probably is not an unmodified face on every"
            << " disconnected part of the patch" << endl
            << "Please use geometric ordering or adapt your topology changes"
            << endl
            << "Continuing with incorrect face ordering from now on!" << endl;
    }
}


bool processorPolyPatch::topologicalOrder
(
    const polyTopoChange& ref,
    const mapPolyMesh& map,
    labelList& faceMap,
    labelList& rotation
) const
{
    const labelList& oldToNew = map.reverseFaceMap();
    const labelList& oldPatchStarts = map.oldPatchStarts();

    const bool isMaster = Pstream::myProcNo() < neighbProcNo();

    // Receive walked faces from neighbour
    labelList masterUnmodOldPatchFaces;
    labelListList masterVisitOrder;
    labelListList masterIndex;
    {
        IPstream fromNeighbour(neighbProcNo());
        fromNeighbour >> masterUnmodOldPatchFaces >> masterVisitOrder
            >> masterIndex;
    }
    if (polyMesh::morphDebug)
    {
        Sout<< "processorPolyPatch::order :"
            << " Received masterSeedFaces:" << masterUnmodOldPatchFaces.size()
            << endl;
    }


    // Get modified/unmodified status for every face on old mesh.
    boolList unmodifiedOldFace(getUnmodifiedFaces(ref, map));

    // Get old patch label from the unmodified faces.
    label oldPatchI = -1;

    forAll(unmodifiedOldFace, oldFaceI)
    {
        if (unmodifiedOldFace[oldFaceI] && inPatch(oldToNew, oldFaceI))
        {
            // Unmodified face in the old patch (since was in
            // processorpatch)
            oldPatchI = whichPatch(oldPatchStarts, oldFaceI);

            if (oldPatchI == -1)
            {
                FatalErrorIn
                (
                    "processorPolyPatch::topologicalOrder"
                    "(const polyTopoChange&, const mapPolyMesh&"
                    ", labelList& faceMap, labelList& rotation)"
                )   << "Unmodified face " << oldToNew[oldFaceI]
                    << " on processor patch " << name()
                    << " was not in a patch in old mesh"
                    << abort(FatalError);
            }

            break;
        }
    }


    // Visited faces.
    boolList visited(size(), false);

    // All patch is one region
    const labelList ppZones(size(), 0);

    bool changed = false;

    forAll(masterUnmodOldPatchFaces, seedI)
    {
        // Patch face (unmodified on master side)
        label masterFaceI = masterUnmodOldPatchFaces[seedI];

        // Old mesh face label
        label oldFaceI = oldPatchStarts[oldPatchI] + masterFaceI;

        if (!unmodifiedOldFace[oldFaceI])
        {
            SeriousError
                << "processorPolyPatch::order on patch:" << name() << " : "
                << "Face " << oldToNew[oldFaceI]
                << " on not modified on processor patch " << name()
                << " but is on its neighbour patch on processor "
                << neighbProcNo() << endl
                << "This means that topological walking on both sides"
                << " will not produce the same ordering of faces" << endl
                << "Please use geometric ordering or adapt your topology"
                << " changes" << endl
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }

        // Walk on patch starting from seed faces to get relative
        // connectivity. Note that seedFaces are the same as on the
        // neighbouring processor since in patch face labels.

        label startFaceI = oldToNew[oldFaceI] - start();

        walkPatch patchWalker
        (
            *this,
            ppZones,
            !isMaster,
            startFaceI,
            localFaces()[startFaceI][0],
            visited
        );

        // Visited path:
        const DynamicList<label>& visitOrder = patchWalker.visitOrder();
        const DynamicList<label>& indexInFace = patchWalker.indexInFace();


        if (polyMesh::morphDebug)
        {
            if (visitOrder.size() > 0 && masterVisitOrder[seedI].size() > 0)
            {
                Sout<< "processorPolyPatch::order : "
                    << "On master walked from unmodified face " << masterFaceI
                    << " and visited " << masterVisitOrder[seedI].size()
                    << " On slave walked from unmodified face " << startFaceI
                    << " and visited " << visitOrder.size() << endl;

                pointField ctrs(calcFaceCentres(*this, map.preMotionPoints()));

                fileName ccName(name() + '_' + Foam::name(seedI) + ".obj");

                Sout<< "processorPolyPatch::order : Dumping visit order for"
                    << " seed " << startFaceI << " to OBJ file " << ccName
                    << endl;

                OFstream ccStr(ccName);
                writeOBJ(ccStr, ctrs, visitOrder);
            }
        }


        // My order should be the same as the one on the master.

        forAll(visitOrder, i)
        {
            label faceI = visitOrder[i];

            label newFaceI = masterVisitOrder[seedI][i];

            faceMap[faceI] = newFaceI;

            // Reorient face
            const face& f = localFaces()[faceI];

            label nShift =
                (
                    (f.size()-masterIndex[seedI][i])
                  - indexInFace[i]
                ) % f.size();

            if (faceI != newFaceI || nShift != 0)
            {
                changed = true;
            }

            rotation[newFaceI] = nShift;
        }
    }

    label falseIndex = findIndex(visited, false);

    if (falseIndex != -1)
    {
        SeriousError
            << "processorPolyPatch::order on patch:" << name() << " : "
            << "Did not visit mesh face " << start() + falseIndex
            << " on processor patch " << name()
            << " from any of the seed faces." << endl
            << "There probably is not an unmodified face on every"
            << " disconnected part of the patch" << endl
            << "Please use geometric ordering or adapt your topology changes"
            << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }

    return changed;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Initialize ordering (on new mesh)
void processorPolyPatch::sendOrder
(
    const polyTopoChange& ref,
    const mapPolyMesh& map
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    const bool isMaster = Pstream::myProcNo() < neighbProcNo();

    if (polyMesh::morphDebug)
    {
        Sout.prefix() = '[' + name() + "] ";
    }

    if (polyMesh::morphDebug)
    {
        // Dump face centres as vertices.
        fileName ctrName(name() + "_ctrs.obj");

        Sout<< "processorPolyPatch::sendOrder : Dumping faceCentres"
            << " to OBJ file " << ctrName
            << " in patch face order" << endl;

        // Recalculate face centres with preMotionPoints (since added faces
        // will have zero area if using old (=current) points)
        pointField ctrs(calcFaceCentres(*this, map.preMotionPoints()));

        OFstream ccStr(ctrName);

        forAll(ctrs, faceI)
        {
            writeOBJ(ccStr, ctrs[faceI]);
        }
    }


    if (geometricMatch() && (separated() || !parallel()))
    {
        SeriousError
            << "processorPolyPatch::order on patch:" << name() << " : "
            << "Cannot use geometric matching on processor patch " << name()
            << " since has 'separated' faces (probably contains part of a"
            << " cyclic patch)" << endl
            << "Please use topological matching instead"
            << " (requires at least on unmodified face per connected area)"
            << " on processor patch " << name()
            << " from any of the seed faces." << endl
            << "Continuing with incorrect face ordering from now on!" << endl;
    }


    if (isMaster)
    {
        if (geometricMatch())
        {
            sendGeometricOrder(ref, map);
        }
        else
        {
            sendTopologicalOrder(ref, map);
        }
    }
}


//- Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool processorPolyPatch::order
(
    const polyTopoChange& ref,
    const mapPolyMesh& map,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    faceMap.setSize(size());
    faceMap = -1;

    rotation.setSize(size());
    rotation = 0;

    const bool isMaster = Pstream::myProcNo() < neighbProcNo();

    bool changed = false;

    if (isMaster)
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFaceI)
        {
            faceMap[patchFaceI] == patchFaceI;
        }
    }
    else
    {
        if (geometricMatch())
        {
            changed = geometricOrder(ref, map, faceMap, rotation);
        }
        else
        {
            changed = topologicalOrder(ref, map, faceMap, rotation);
        }
    }

    return changed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
