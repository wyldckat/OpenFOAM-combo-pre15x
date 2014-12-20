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

\*---------------------------------------------------------------------------*/

#include "coupledPolyPatch.H"
#include "SortableList.H"
#include "ListOps.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledPolyPatch, 0);

scalar coupledPolyPatch::matchTol_ = 1E-3;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void coupledPolyPatch::writeOBJ(Ostream& os, const point& pt)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
}


void coupledPolyPatch::writeOBJ
(
    Ostream& os,
    const pointField& points,
    const labelList& pointLabels
)
{
    forAll(pointLabels, i)
    {
        writeOBJ(os, points[pointLabels[i]]);
    }
}


// Write edge
void coupledPolyPatch::writeOBJ
(
    Ostream& os,
    const point& p0,
    const point& p1,
    label& vertI
)
{
    writeOBJ(os, p0);
    vertI++;

    writeOBJ(os, p1);
    vertI++;

    os<< "l " << vertI-1 << ' ' << vertI << nl;
}


// Get face centre with specified points
// (can't use points() since might not be set yet in morphing; usually called
//  with mapPolyMesh.premotionPoints())
pointField coupledPolyPatch::calcFaceCentres
(
    const faceList& faces,
    const pointField& points
)
{
    pointField ctrs(faces.size());

    forAll(faces, faceI)
    {
        ctrs[faceI] = faces[faceI].centre(points);
    }

    return ctrs;
}


// Get coordinate of f[0] for every face
pointField coupledPolyPatch::getAnchorPoints
(
    const faceList& faces,
    const pointField& points
)
{
    pointField anchors(faces.size());

    forAll(faces, faceI)
    {
        anchors[faceI] = points[faces[faceI][0]];
    }

    return anchors;
}


// Is old face now in the current patch?
bool coupledPolyPatch::inPatch
(
    const labelList& oldToNew,
    const label oldFaceI
) const
{
    label faceI = oldToNew[oldFaceI];

    return faceI >= start() && faceI < start()+size();
}


// Given list of starts of patches and a face label determine the patch.
label coupledPolyPatch::whichPatch
(
    const labelList& patchStarts,
    const label faceI
)
{
    forAll(patchStarts, patchI)
    {
        if (patchStarts[patchI] <= faceI)
        {
            if (patchI == patchStarts.size()-1)
            {
                return patchI;
            }
            else if (patchStarts[patchI+1] > faceI)
            {
                return patchI;
            }
        }
    }

    return -1;
}


// Get local typical dimension and tolerance from that. Currently max of
// distance from centre to any of the face points.
scalarField coupledPolyPatch::calcFaceTol
(
    const faceList& faces,
    const pointField& points,
    const pointField& faceCentres
)
{
    // Calculate typical distance per face
    scalarField tols(faces.size());

    forAll(faces, faceI)
    {
        const point& cc = faceCentres[faceI];

        const face& f = faces[faceI];

        scalar maxLen = -GREAT;

        forAll(f, fp)
        {
            maxLen = max(maxLen, mag(points[f[fp]] - cc));
        }
        tols[faceI] = matchTol_ * maxLen;
    }
    return tols;
}


label coupledPolyPatch::getRotation
(
    const pointField& points,
    const face& f,
    const point& anchor,
    const scalar tol
)
{
    label anchorFp = -1;
    scalar minDistSqr = GREAT;

    forAll(f, fp)
    {
        scalar distSqr = magSqr(anchor - points[f[fp]]);

        if (distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            anchorFp = fp;
        }
    }

    if (anchorFp == -1 || mag(minDistSqr) > tol)
    {
        return -1;
    }
    else
    {
        // Positive rotation
        return (f.size() - anchorFp) % f.size();
    }
}


void coupledPolyPatch::calcTransformTensors
(
    const vector& Cf,
    const vector& Cr,
    const vector& nf,
    const vector& nr
) const
{
    if (mag(nf & nr) < 1 - SMALL)
    {
        separation_.setSize(0);

        forwardT_ = tensorField(1, transformationTensor(-nr, nf));
        reverseT_ = tensorField(1, transformationTensor(nf, -nr));
    }
    else
    {
        forwardT_.setSize(0);
        reverseT_.setSize(0);

        vector separation = (nf & (Cr - Cf))*nf;

        if (mag(separation) > SMALL)
        {
            separation_ = vectorField(1, separation);
        }
        else
        {
            separation_.setSize(0);
        }
    }
}


void coupledPolyPatch::calcTransformTensors
(
    const vectorField& Cf,
    const vectorField& Cr,
    const vectorField& nf,
    const vectorField& nr
) const
{
    if (size() == 0)
    {
        // Dummy geometry.
        separation_.setSize(0);
        forwardT_ = I;
        reverseT_ = I;
    }
    else if (sum(mag(nf & nr)) < size() - SMALL)
    {
        separation_.setSize(0);

        forwardT_.setSize(size());
        reverseT_.setSize(size());

        forAll (forwardT_, facei)
        {
            forwardT_[facei] = transformationTensor(-nr[facei], nf[facei]);
            reverseT_[facei] = transformationTensor(nf[facei], -nr[facei]);
        }

        if (sum(mag(forwardT_ - forwardT_[0])) < SMALL)
        {
            forwardT_.setSize(1);
            reverseT_.setSize(1);
        }
    }
    else
    {
        forwardT_.setSize(0);
        reverseT_.setSize(0);

        separation_ = (nf&(Cr - Cf))*nf;

        if (sum(mag(separation_))/size() < SMALL)
        {
            separation_.setSize(0);
        }
        else if (sum(mag(separation_ - separation_[0]))/size() < SMALL)
        {
            separation_.setSize(1);
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm)
{}


coupledPolyPatch::coupledPolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(is, index, bm)
{}


coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm)
{}


coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm)
{}


coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coupledPolyPatch::~coupledPolyPatch()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
