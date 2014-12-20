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

\*---------------------------------------------------------------------------*/

#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicPolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, word);
addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, Istream);
addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Force calculation of transformation tensors
void cyclicPolyPatch::calcTransforms()
{
    if (size() > 0)
    {
        const pointField& points = allPoints();

        const face& f0 = static_cast<const faceList&>(*this)[0];
        const face& fn2 = static_cast<const faceList&>(*this)[size()/2];

        vector nf0 = f0.normal(points);
        nf0 /= mag(nf0);

        vector nfn2 = fn2.normal(points);
        nfn2 /= mag(nfn2);

        calcTransformTensors
        (
            f0.centre(points),
            fn2.centre(points),
            nf0,
            nfn2
        );
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

// Construct from components
cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}


// Construct from Istream
cyclicPolyPatch::cyclicPolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(is, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}


// Construct from dictionary
cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}

//- Construct as copy, resetting the boundary mesh
cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}

//- Construct as copy, resetting the face list and boundary mesh data
cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cyclicPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}

void cyclicPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}

void cyclicPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}

void cyclicPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    calcTransforms();
}

void cyclicPolyPatch::initUpdateTopology()
{
    polyPatch::initUpdateTopology();
}

void cyclicPolyPatch::updateTopology()
{
    polyPatch::updateTopology();
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


const edgeList& cyclicPolyPatch::coupledPoints() const
{
    if (!coupledPointsPtr_)
    {
        // Look at cyclic patch as two halves, A and B.
        // Now all we know is that relative face index in halfA is same
        // as coupled face in halfB and also that the 0th vertex
        // corresponds.

        // From halfA point to halfB or -1.
        labelList coupledPoint(nPoints(), -1);

        for (label patchFaceA = 0; patchFaceA < size()/2; patchFaceA++)
        {
            const face& fA = localFaces()[patchFaceA];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (coupledPoint[patchPointA] == -1)
                {
                    const face& fB = localFaces()[patchFaceA + size()/2];

                    label indexB = (fB.size() - indexA) % fB.size();

                    coupledPoint[patchPointA] = fB[indexB];
                }
            }
        }

        coupledPointsPtr_ = new edgeList(nPoints());
        edgeList& connected = *coupledPointsPtr_;

        // Extract coupled points.
        label connectedI = 0;

        forAll(coupledPoint, i)
        {
            if (coupledPoint[i] != -1)
            {
                connected[connectedI++] = edge(i, coupledPoint[i]);
            }
        }

        connected.setSize(connectedI);

        if (debug)
        {
            Info<< "Writing file coupledPoints.obj with coordinates of "
                << "coupled points" << endl;

            OFstream str("coupledPoints.obj");
            label vertI = 0;

            forAll(connected, i)
            {
                const point& a = localPoints()[connected[i][0]];
                const point& b = localPoints()[connected[i][1]];

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledPointsPtr_;
}


const edgeList& cyclicPolyPatch::coupledEdges() const
{
    if (!coupledEdgesPtr_)
    {
        // Build map from points on halfA to points on halfB.
        const edgeList& pointCouples = coupledPoints();

        Map<label> aToB(2*pointCouples.size());

        forAll(pointCouples, i)
        {
            const edge& e = pointCouples[i];

            aToB.insert(e[0], e[1]);
        }

        // Map from edge on half A to points (in halfB indices)
        HashTable<label, edge, Hash<edge> > edgeMap(nEdges());

        for (label patchFaceA = 0; patchFaceA < size()/2; patchFaceA++)
        {
            const labelList& fEdges = faceEdges()[patchFaceA];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Convert edge end points to corresponding points on halfB
                // side.

                edge halfBEdge(aToB[e[0]], aToB[e[1]]);

                edgeMap.insert(halfBEdge, edgeI);
            }
        }

        coupledEdgesPtr_ = new edgeList(nEdges()/2);
        edgeList& coupledEdges = *coupledEdgesPtr_;
        label coupleI = 0;

        for (label patchFaceB = size()/2; patchFaceB < size(); patchFaceB++)
        {
            const labelList& fEdges = faceEdges()[patchFaceB];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Look up edge from HashTable.
                label halfAEdgeI = edgeMap[e];

                // Store correspondence
                coupledEdges[coupleI++] = edge(halfAEdgeI, edgeI);

                // Remove (since should be used only once)
                edgeMap.erase(e);
            }
        }


        // Some checks

        if (coupleI != coupledEdges.size())
        {
            FatalErrorIn("cyclicPolyPatch::coupledEdges() const")
                << "Problem : coupleI:" << coupleI
                << " number of edges per halfpatch:" << coupledEdges.size()
                << abort(FatalError);
        }

        forAll(coupledEdges, i)
        {
            const edge& e = coupledEdges[i];

            if (e[0] == e[1] || e[0] < 0 || e[1] < 0)
            {
                FatalErrorIn("cyclicPolyPatch::coupledEdges() const")
                    << "Problem : at position " << i
                    << " illegal couple:" << e
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            Info<< "Writing file coupledEdges.obj with centres of "
                << "coupled edges" << endl;

            OFstream str("coupledEdges.obj");
            label vertI = 0;

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                const point& a = edges()[e[0]].centre(localPoints());
                const point& b = edges()[e[1]].centre(localPoints());

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledEdgesPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
