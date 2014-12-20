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

#include "error.H"

#include "pointBoundaryMesh.H"
#include "polyBoundaryMesh.H"
#include "facePointPatch.H"
#include "globalProcessorPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyBoundaryMesh
pointBoundaryMesh::pointBoundaryMesh
(
    const pointMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    pointPatchList(basicBdry.size()),
    mesh_(m)
{
    // hook boundary patches
    pointPatchList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.hook(facePointPatch::New(basicBdry[patchI], *this).ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointPatch& pointBoundaryMesh::globalPointPatch() const
{
    const pointPatchList& patches = *this;

    forAll (patches, patchI)
    {
        if
        (
            typeid(patches[patchI])
         == typeid(globalProcessorPointPatch)
        )
        {
            return patches[patchI];
        }
    }

    FatalErrorIn
    (
        "const pointBoundaryMesh::"
        "globalProcessorPointPatch& globalPointPatch() const"
    )   << "patch not found.  Is this case running in parallel?"
        << abort(FatalError);

    // Dummy return
    return patches[0];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
