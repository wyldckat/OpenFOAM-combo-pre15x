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

Class
    syncTools

Description

\*----------------------------------------------------------------------------*/

#include "syncTools.H"
#include "PstreamReduceOps.H"
#include "PstreamCombineReduceOps.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "parallelInfo.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Foam::syncTools::transformList
(
    const tensorField& rotTensor,
    UList<T>& field
)
{
    if (rotTensor.size() == 1)
    {
        forAll(field, i)
        {
            field[i] = Foam::transform(rotTensor[0], field[i]);
        }
    }
    else if (rotTensor.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] = Foam::transform(rotTensor[i], field[i]);
        }
    }
    else
    {
        FatalErrorIn
        (
            "syncTools::separateList(const vectorField&, UList<T>&)"
        )   << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


template <class T>
void Foam::syncTools::separateList
(
    const vectorField& separation,
    UList<T>& field
)
{}


template <class T, class CombineOp>
Foam::label Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    UList<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const bool applySeparation
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncPointList"
            "(const polyMesh&, UList<T>&, const CombineOp&, const T&"
            ", const bool)"
        )   << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return 0;
    }

    label nChanged = 0;

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            const polyPatch& patch = patches[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Get data per patchPoint in neighbouring point numbers.
                List<T> patchInfo(procPatch.nPoints());

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll(meshPts, pointI)
                {
                    patchInfo[nbrPts[pointI]] = pointValues[meshPts[pointI]];
                }

                OPstream toNeighb(procPatch.neighbProcNo());
                toNeighb << patchInfo;
            }
        }

        // Receive and combine.

        forAll(patches, patchI)
        {
            const polyPatch& patch = patches[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                checkTransform(procPatch, applySeparation);

                List<T> nbrPatchInfo(procPatch.nPoints());
                {
                    IPstream fromNeighb(procPatch.neighbProcNo());
                    fromNeighb >> nbrPatchInfo;
                }

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.reverseT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointI)
                {
                    label meshPointI = meshPts[pointI];

                    const T oldVal(pointValues[meshPointI]);

                    cop(pointValues[meshPointI], nbrPatchInfo[pointI]);

                    if (oldVal != pointValues[meshPointI])
                    {
                        nChanged++;
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            checkTransform(cycPatch, applySeparation);

            const labelList& meshPts = cycPatch.meshPoints();
            const edgeList& coupledPoints = cycPatch.coupledPoints();

            List<T> half0Values(coupledPoints.size());
            List<T> half1Values(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                label point0 = meshPts[e[0]];
                label point1 = meshPts[e[1]];

                half0Values[i] = pointValues[point0];
                half1Values[i] = pointValues[point1];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.forwardT(), half0Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                hasTransformation = true;

                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                label point0 = meshPts[e[0]];
                label point1 = meshPts[e[1]];

                const T oldVal(pointValues[point0]);

                cop(pointValues[point0], half1Values[i]);
                cop(pointValues[point1], half0Values[i]);

                if (oldVal != pointValues[point0])
                {
                    nChanged++;
                }
            }
        }
    }


    // Synchronize multiple shared points.
    const parallelInfo& pd = mesh.parallelData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningIn
            (
                "syncTools<class T, class CombineOp>::syncPointList"
                "(const polyMesh&, UList<T>&, const CombineOp&, const T&"
                ", const bool)"
            )   << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }
        else
        {
            // Values on shared points.
            List<T> sharedPts(pd.nGlobalPoints(), nullValue);

            forAll(pd.sharedPointLabels(), i)
            {
                label meshPointI = pd.sharedPointLabels()[i];

                // Fill my entries in the shared points
                sharedPts[pd.sharedPointAddr()[i]] = pointValues[meshPointI];
            }

            // Combine on master.
            Pstream::listCombineGather(sharedPts, cop);
            Pstream::listCombineScatter(sharedPts);

            // Now we will all have the same information. Merge it back with
            // my local information.
            forAll(pd.sharedPointLabels(), i)
            {
                label meshPointI = pd.sharedPointLabels()[i];

                const T& sharedVal = sharedPts[pd.sharedPointAddr()[i]];

                if (sharedVal != pointValues[meshPointI])
                {
                    pointValues[meshPointI] = sharedVal;
                    nChanged++;
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());

    return nChanged;
}


template <class T, class CombineOp>
Foam::label Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    UList<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const bool applySeparation
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncEdgeList"
            "(const polyMesh&, UList<T>&, const CombineOp&, const T&"
            ", const bool)"
        )   << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return 0;
    }

    label nChanged = 0;

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Collect and send information.
        // Send region info for all patch edges
        forAll(patches, patchI)
        {
            if (isA<processorPolyPatch>(patches[patchI]))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Get region per patch edge in neighbouring edge numbers.
                List<T> patchInfo(procPatch.nEdges());

                const labelList& meshEdges = procPatch.meshEdges();
                const labelList& neighbEdges = procPatch.neighbEdges();

                forAll(meshEdges, edgeI)
                {
                    patchInfo[neighbEdges[edgeI]] =
                        edgeValues[meshEdges[edgeI]];
                }

                OPstream toNeighb(procPatch.neighbProcNo());
                toNeighb << patchInfo;
            }
        }

        // Receive and merge.

        forAll(patches, patchI)
        {
            if (isA<processorPolyPatch>(patches[patchI]))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                checkTransform(procPatch, applySeparation);

                const labelList& meshEdges = procPatch.meshEdges();

                // Receive from neighbour. Is per patch edge the region of the
                // neighbouring patch edge. 
                List<T> nbrPatchInfo(procPatch.nEdges());
                {
                    IPstream fromNeighb(procPatch.neighbProcNo());
                    fromNeighb >> nbrPatchInfo;
                }

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.reverseT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }


                forAll(meshEdges, edgeI)
                {
                    label meshEdgeI = meshEdges[edgeI];

                    const T oldVal(edgeValues[meshEdgeI]);

                    cop(edgeValues[meshEdgeI], nbrPatchInfo[edgeI]);

                    if (oldVal != edgeValues[meshEdgeI])
                    {
                        nChanged++;
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            checkTransform(cycPatch, applySeparation);

            const edgeList& coupledEdges = cycPatch.coupledEdges();
            const labelList& meshEdges = cycPatch.meshEdges();

            List<T> half0Values(coupledEdges.size());
            List<T> half1Values(coupledEdges.size());

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                label meshEdge0 = meshEdges[e[0]];
                label meshEdge1 = meshEdges[e[1]];

                half0Values[i] = edgeValues[meshEdge0];
                half1Values[i] = edgeValues[meshEdge1];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.forwardT(), half0Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                hasTransformation = true;

                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                label meshEdge0 = meshEdges[e[0]];
                label meshEdge1 = meshEdges[e[1]];

                const T oldVal(edgeValues[meshEdge0]);

                cop(edgeValues[meshEdge0], half1Values[i]);
                cop(edgeValues[meshEdge1], half0Values[i]);

                if (oldVal != edgeValues[meshEdge0])
                {
                    nChanged++;
                }
            }
        }
    }

    // Do the multiple shared edges
    const parallelInfo& pd = mesh.parallelData();

    if (pd.nGlobalEdges() > 0)
    {
        if (hasTransformation)
        {
            WarningIn
            (
                "syncTools<class T, class CombineOp>::syncEdgeList"
                "(const polyMesh&, UList<T>&, const CombineOp&, const T&"
                ", const bool)"
            )   << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }
        else
        {
            // Values on shared edges.
            List<T> sharedPts(pd.nGlobalEdges(), nullValue);

            forAll(pd.sharedEdgeLabels(), i)
            {
                label meshEdgeI = pd.sharedEdgeLabels()[i];

                // Fill my entries in the shared edges
                sharedPts[pd.sharedEdgeAddr()[i]] = edgeValues[meshEdgeI];
            }

            // Combine on master.
            Pstream::listCombineGather(sharedPts, cop);
            Pstream::listCombineScatter(sharedPts);

            // Now we will all have the same information. Merge it back with
            // my local information.
            forAll(pd.sharedEdgeLabels(), i)
            {
                label meshEdgeI = pd.sharedEdgeLabels()[i];

                const T& sharedVal = sharedPts[pd.sharedEdgeAddr()[i]];

                if (sharedVal != edgeValues[meshEdgeI])
                {
                    edgeValues[meshEdgeI] = sharedVal;
                    nChanged++;
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());

    return nChanged;
}


template <class T, class CombineOp>
Foam::label Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const bool applySeparation
)
{
    if (faceValues.size() != mesh.nFaces())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncFaceList"
            "(const polyMesh&, UList<T>&, const CombineOp&"
            ", const bool)"
        )   << "Number of values " << faceValues.size()
            << " is not equal to the number of faces in the mesh "
            << mesh.nFaces() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return 0;
    }


    label nChanged = 0;

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            const polyPatch& patch = patches[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                OPstream toNeighb(procPatch.neighbProcNo());
                toNeighb <<
                    SubList<T>(faceValues, procPatch.size(), procPatch.start());
            }
        }

        // Receive and combine.

        forAll(patches, patchI)
        {
            const polyPatch& patch = patches[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<T> nbrPatchInfo(procPatch.size());
                {
                    IPstream fromNeighb(procPatch.neighbProcNo());
                    fromNeighb >> nbrPatchInfo;
                }

                if (!procPatch.parallel())
                {
                    transformList(procPatch.reverseT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    separateList(procPatch.separation(), nbrPatchInfo);
                }

                label meshFaceI = procPatch.start();

                forAll(nbrPatchInfo, faceI)
                {
                    const T oldVal(faceValues[meshFaceI]);

                    cop(faceValues[meshFaceI], nbrPatchInfo[faceI]);

                    if (oldVal != faceValues[meshFaceI])
                    {
                        nChanged++;
                    }
                    meshFaceI++;
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            label half = cycPatch.size()/2;
            label half1Start = cycPatch.start()+half;

            List<T> half0Values(SubList<T>(faceValues, half, cycPatch.start()));
            List<T> half1Values(SubList<T>(faceValues, half, half1Start));

            // Save old values.
            const List<T> oldVal(half0Values);

            if (!cycPatch.parallel())
            {
                transformList(cycPatch.reverseT(), half1Values);
                transformList(cycPatch.forwardT(), half0Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            label mesh0FaceI = cycPatch.start();
            label mesh1FaceI = cycPatch.start() + half;

            forAll(half0Values, faceI)
            {
                cop(faceValues[mesh0FaceI++], half1Values[faceI]);
                cop(faceValues[mesh1FaceI++], half0Values[faceI]);
            }

            mesh0FaceI = cycPatch.start();

            forAll(oldVal, faceI)
            {
                if (faceValues[mesh0FaceI++] != oldVal[faceI])
                {
                    nChanged++;
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());

    return nChanged;
}


template <class T>
void Foam::syncTools::swapFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const bool applySeparation
)
{
    if (faceValues.size() != mesh.nFaces())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::swapFaceList"
            "(const polyMesh&, UList<T>&, const bool)"
        )   << "Number of values " << faceValues.size()
            << " is not equal to the number of faces in the mesh "
            << mesh.nFaces() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }


    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            const polyPatch& patch = patches[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                OPstream toNeighb(procPatch.neighbProcNo());
                toNeighb <<
                    SubList<T>(faceValues, procPatch.size(), procPatch.start());
            }
        }

        // Receive.

        forAll(patches, patchI)
        {
            const polyPatch& patch = patches[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<T> nbrPatchInfo(procPatch.size());
                {
                    IPstream fromNeighb(procPatch.neighbProcNo());
                    fromNeighb >> nbrPatchInfo;
                }

                if (!procPatch.parallel())
                {
                    transformList(procPatch.reverseT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    separateList(procPatch.separation(), nbrPatchInfo);
                }


                label meshFaceI = procPatch.start();

                forAll(nbrPatchInfo, i)
                {
                    faceValues[meshFaceI++] = nbrPatchInfo[i];
                }

            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            label half = cycPatch.size()/2;

            SubList<T> half0(faceValues, half, cycPatch.start());
            SubList<T> half1(faceValues, half, cycPatch.start()+half);

            List<T> half0Values(half0);
            List<T> half1Values(half1);

            if (!cycPatch.parallel())
            {
                transformList(cycPatch.forwardT(), half0Values);
                transformList(cycPatch.reverseT(), half1Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            // Swap values
            label mesh0FaceI = cycPatch.start();
            label mesh1FaceI = cycPatch.start()+half;

            forAll(half0Values, i)
            {
                faceValues[mesh0FaceI++] = half1Values[i];
                faceValues[mesh1FaceI++] = half0Values[i];
            }
        }
    }
}


// ************************************************************************* //
