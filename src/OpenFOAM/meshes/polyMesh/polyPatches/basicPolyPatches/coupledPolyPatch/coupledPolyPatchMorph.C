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

#include "coupledPolyPatch.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "SortableList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool coupledPolyPatch::geometricMatch_ = true;

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


// Get unmodified faces on original mesh.
boolList coupledPolyPatch::getUnmodifiedFaces
(
    const polyTopoChange& ref,
    const mapPolyMesh& map
)
{
    boolList unmodifiedFace(map.nOldFaces(), true);

    // Unmark any modified faces
    const DynamicList<polyModifyFace>& mf = ref.modifiedFaces();

    forAll (mf, mfI)
    {
        label oldFaceI = mf[mfI].faceID();

        unmodifiedFace[oldFaceI] = false;
    }

    // Any removed faces
    const labelHashSet& removedFaces = ref.removedFaces();

    for
    (
        labelHashSet::const_iterator iter = removedFaces.begin();
        iter != removedFaces.end();
        ++iter
    )
    {
        unmodifiedFace[iter.key()] = false;
    }
    return unmodifiedFace;
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
) const
{
    label anchorFp = -1;

    forAll(f, fp)
    {
        if (mag(anchor - points[f[fp]]) < tol)
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledPolyPatch::setGeometricMatch(const bool geometricMatch)
{
    geometricMatch_ = geometricMatch;
}


bool coupledPolyPatch::geometricMatch()
{
    return geometricMatch_;
}


void coupledPolyPatch::setGeometricMatchTol(const scalar matchTol)
{
    matchTol_ = matchTol;
}


scalar coupledPolyPatch::geometricMatchTol()
{
    return matchTol_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
