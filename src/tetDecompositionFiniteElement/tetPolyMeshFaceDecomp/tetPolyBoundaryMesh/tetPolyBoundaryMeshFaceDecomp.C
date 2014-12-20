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

#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "polyBoundaryMesh.H"
#include "faceTetPolyPatchFaceDecomp.H"
#include "globalProcessorTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyBoundaryMesh
tetPolyBoundaryMeshFaceDecomp::tetPolyBoundaryMeshFaceDecomp
(
    const tetPolyMeshFaceDecomp& m,
    const polyBoundaryMesh& basicBdry
)
:
    tetPolyPatchFaceDecompList(basicBdry.size()),
    mesh_(m)
{
    // Hook boundary patches
    tetPolyPatchFaceDecompList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.hook
        (
            faceTetPolyPatchFaceDecomp::New(basicBdry[patchI], *this).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const tetPolyPatchFaceDecomp&
tetPolyBoundaryMeshFaceDecomp::globalPointPatch() const
{
    const tetPolyPatchFaceDecompList& patches = *this;

    forAll (patches, patchI)
    {
        if
        (
            typeid(patches[patchI])
         == typeid(globalProcessorTetPolyPatchFaceDecomp)
        )
        {
            return patches[patchI];
        }
    }

    FatalErrorIn
    (
        "const tetPolyBoundaryMeshFaceDecomp::"
        "globalProcessorTetPolyPatchFaceDecomp& globalPointPatch() const"
    )   << "patch not found.  Is this case running in parallel?"
        << abort(FatalError);

    // Dummy return
    return patches[0];
}


faceListList tetPolyBoundaryMeshFaceDecomp::boundaryTriFaces() const
{
    faceListList result(size());

    forAll (result, patchI)
    {
        result[patchI] = operator[](patchI).triFaces();
    }

    return result;
}


void tetPolyBoundaryMeshFaceDecomp::updateMesh()
{
    // Hook boundary patches
    tetPolyPatchFaceDecompList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches[patchI].updateMesh();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
