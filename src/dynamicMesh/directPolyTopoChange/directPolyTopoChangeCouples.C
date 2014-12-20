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
    Almost straight copy of polyMeshMorph using geometric matching.

\*---------------------------------------------------------------------------*/

#include "directPolyTopoChange.H"
#include "Time.H"
#include "morphMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get face centre for slice of faces.
Foam::pointField Foam::directPolyTopoChange::calcFaceCentres
(
    const label start,
    const label size
) const
{
    pointField ctrs(size, vector::zero);

    forAll(ctrs, i)
    {
        label faceI = start + i;

        const face& f = faces_[faceI];

        forAll(f, fp)
        {
            ctrs[i] += points_[f[fp]];
        }
        ctrs[i] /= f.size();
    }

    return ctrs;
}


// Get coordinate of f[0] for slice of faces
Foam::pointField Foam::directPolyTopoChange::getAnchorPoints
(
    const label start,
    const label size
) const
{
    pointField anchors(size);

    forAll(anchors, i)
    {
        label faceI = start + i;

        anchors[i] = points_[faces_[faceI][0]];
    }

    return anchors;
}


// Get local typical dimension and tolerance from that. Currently max of
// distance from centre to any of the face points.
Foam::scalarField Foam::directPolyTopoChange::calcFaceTol
(
    const label start,
    const label size,
    const pointField& faceCentres
) const
{
    // Calculate typical distance per face
    scalarField tols(size);

    forAll(tols, i)
    {
        label faceI = i + start;

        const point& cc = faceCentres[i];

        const face& f = faces_[faceI];

        scalar maxLen = -GREAT;

        forAll(f, fp)
        {
            maxLen = max(maxLen, mag(points_[f[fp]] - cc));
        }
        tols[i] = 1E-3 * maxLen;
    }
    return tols;
}

        
Foam::label Foam::directPolyTopoChange::getRotation
(
    const face& f,
    const point& anchor,
    const scalar tol
) const
{
    label anchorFp = -1;

    forAll(f, fp)
    {
        if (mag(anchor - points_[f[fp]]) < tol)
        {
            anchorFp = fp;
            break;
        }
    }

    if (anchorFp == -1)
    {
        return -1;
    }
    else
    {
        // Positive rotation
        return (f.size() - anchorFp) % f.size();
    }
}


Foam::face Foam::directPolyTopoChange::rotateFace
(
    const face& f,
    const label nPos
)
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


void Foam::directPolyTopoChange::initOrderProcessorPolyPatch
(
    const label neighbProcNo,
    const label start,
    const label size
) const
{
    if (Pstream::myProcNo() < neighbProcNo)
    {
        Pout<< "Sending to " << neighbProcNo<< endl
            << "    faceCentres:" << calcFaceCentres(start, size) << endl
            << "    anchors    :" << getAnchorPoints(start, size) << endl;

        OPstream toNeighbour(neighbProcNo);
        toNeighbour
            << calcFaceCentres(start, size)
            << getAnchorPoints(start, size);
    }
}


void Foam::directPolyTopoChange::initOrderCyclicPolyPatch
(
    const label start,
    const label size
) const
{
    SeriousError << "directPolyTopoChange::initOrderCyclicPolyPatch : "
        << "Cannot handle cyclics" << endl;
}


bool Foam::directPolyTopoChange::orderProcessorPolyPatch
(
    const label neighbProcNo,
    const label start,
    const label size,

    labelList& faceMap,
    labelList& rotation
) const
{
    if (Pstream::myProcNo() < neighbProcNo)
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        faceMap = identity(faceMap.size());

        return false;   // nothing changed
    }
    else
    {
        pointField masterCtrs;
        pointField masterAnchorPoints;

        // Receive data from neighbour
        {
            IPstream fromNeighbour(neighbProcNo);
            fromNeighbour >> masterCtrs >> masterAnchorPoints;
        }

        Pout<< "Received from " << neighbProcNo<< endl
            << "    faceCentres:" << calcFaceCentres(start, size) << endl
            << "    anchors    :" << getAnchorPoints(start, size) << endl;


        // Calculate my face centres
        pointField ctrs(calcFaceCentres(start, size));

        // Calculate typical distance from face centre
        scalarField tols(calcFaceTol(start, size, ctrs));

        // Geometric match of face centre vectors
        bool matchedAll = matchPoints(ctrs, masterCtrs, tols, false, faceMap);

        if (!matchedAll)
        {
            SeriousError
                << "directPolyTopoChange::orderProcessorPolyPatch :"
                << " Cannot match vectors to faces on both sides of patch"
                << endl
                << "masterCtrs[0]:" << masterCtrs[0] << endl
                << "ctrs[0]:" << ctrs[0] << endl
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
                    faces_[oldFaceI + start],
                    wantedAnchor,
                    tols[oldFaceI]
                );

            if (rotation[newFaceI] == -1)
            {
                SeriousError
                    << "Cannot find point on face " << faces_[oldFaceI + start]
                    << " with vertices:"
                    << IndirectList<point>(points_, faces_[oldFaceI + start])
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch" << endl
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;

                return false;
            }
        }

        Pout<< "faceMap :" << faceMap << endl;
        Pout<< "rotation:" << rotation << endl;

        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] != faceI || rotation[faceI] != 0)
            {
                return true;
            }
        }
        return false;
    }
}


bool Foam::directPolyTopoChange::orderCyclicPolyPatch
(
    const label start,
    const label size,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Do nothing (i.e. identical mapping, zero rotation).
    faceMap = identity(faceMap.size());

    return false;   // nothing changed
}


void Foam::directPolyTopoChange::reorderCoupledFaces
(
    const PtrList<dictionary>& oldPatchDicts,
    const labelList& patchStarts,
    const labelList& patchSizes
)
{
    // Extract patch types
    wordList oldPatchTypes(oldPatchDicts.size());

    forAll(oldPatchDicts, patchI)
    {
        oldPatchTypes[patchI] = word(oldPatchDicts[patchI].lookup("type"));
    }


    // Mapping for faces (old to new). Extends over all mesh faces for
    // convenience (could be just the external faces)
    labelList faceMap(faces_.size());
    forAll(faceMap, faceI)
    {
        faceMap[faceI] = faceI;
    }

    // Rotation on new faces.
    labelList rotation(faces_.size(), 0);


    // Send ordering
    forAll(oldPatchTypes, patchI)
    {
        if 
        (
            Pstream::parRun()
         && oldPatchTypes[patchI] == processorPolyPatch::typeName
        )
        {
            initOrderProcessorPolyPatch
            (
                readLabel(oldPatchDicts[patchI].lookup("neighbProcNo")),
                patchStarts[patchI],
                patchSizes[patchI]
            );
        }
        else if (oldPatchTypes[patchI] == cyclicPolyPatch::typeName)
        {
            initOrderCyclicPolyPatch
            (
                patchStarts[patchI],
                patchSizes[patchI]
            );
        }
    }

    // Receive and calculate ordering

    bool anyChanged = false;

    forAll(oldPatchTypes, patchI)
    {
        labelList patchFaceMap(patchSizes[patchI], -1);
        labelList patchFaceRotation(patchSizes[patchI], 0);

        bool changed = false;

        if 
        (
            Pstream::parRun()
         && oldPatchTypes[patchI] == processorPolyPatch::typeName
        )
        {
            changed = orderProcessorPolyPatch
            (
                readLabel(oldPatchDicts[patchI].lookup("neighbProcNo")),
                patchStarts[patchI],
                patchSizes[patchI],
                patchFaceMap,
                patchFaceRotation
            );
        }
        else if (oldPatchTypes[patchI] == cyclicPolyPatch::typeName)
        {
            changed = orderCyclicPolyPatch
            (
                patchStarts[patchI],
                patchSizes[patchI],
                patchFaceMap,
                patchFaceRotation
            );
        }

        if (changed)
        {
            // Merge patch face reordering into mesh face reordering table
            label start = patchStarts[patchI];

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
        // Reorder faces according to faceMap.
        reorderCompactFaces(faceMap.size(), faceMap, false);

        // Rotate faces (rotation is already in new face indices).
        forAll(rotation, faceI)
        {
            if (rotation[faceI] != 0)
            {
                faces_[faceI] = rotateFace(faces_[faceI], rotation[faceI]);
            }
        }
    }
}


// ************************************************************************* //
