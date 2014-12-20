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

        const face& f0 = ((const faceList&)(*this))[0];
        const face& fn2 = ((const faceList&)(*this))[size()/2];

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
    coupledPointsPtr_(NULL)
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
    coupledPointsPtr_(NULL)
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
    coupledPointsPtr_(NULL)
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
    coupledPolyPatch(pp, bm)
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
    coupledPointsPtr_(NULL)
{
    calcTransforms();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cyclicPolyPatch::initGeometry()
{}

void cyclicPolyPatch::calcGeometry()
{}

void cyclicPolyPatch::initMovePoints(const pointField&)
{}

void cyclicPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    calcTransforms();
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
