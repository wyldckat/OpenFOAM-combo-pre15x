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

#include "facePointPatch.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#include "demandDrivenData.H"
#include "boolList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(facePointPatch, 0);
defineRunTimeSelectionTable(facePointPatch, polyPatch);

addToRunTimeSelectionTable
(
    facePointPatch,
    facePointPatch,
    polyPatch
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyPatch
facePointPatch::facePointPatch
(
    const polyPatch& p,
    const pointBoundaryMesh& bm
)
:
    pointPatch(bm),
    polyPatch_(p)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

triFaceList facePointPatch::faceTriangles
(
    const label faceID
) const
{
    const face& f =
        boundaryMesh().mesh()()
            .boundaryMesh()[index()].localFaces()[faceID];

    // Create a list of triangles to keep the triangles that
    // have already been added
    triFaceList result(f.size() - 2);

    for (label triI = 0; triI < (f.size() - 2); triI++)
    {
        result[triI] = triFace(f[0], f[triI + 1], f[triI + 2]);
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
