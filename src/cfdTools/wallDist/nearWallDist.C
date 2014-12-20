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

#include "nearWallDist.H"
#include "fvMesh.H"
#include "cellDistFuncs.H"
#include "wallFvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDist::doAll()
{
    cellDistFuncs wallUtils(mesh_);

    // Get patch ids of walls
    labelHashSet wallPatchIDs
    (
        wallUtils.getPatchIDs
        (
            wallPolyPatch::typeName
        )
    );

    // Size neighbours array for maximum possible

    labelList neighbours(wallUtils.maxPatchSize(wallPatchIDs));


    // Correct all cells with face on wall

    const volVectorField& cellCentres = mesh_.C();

    forAll(mesh_.boundary(), patchI)
    {
        fvPatchScalarField& ypatch = operator[](patchI);

        const fvPatch& patch = mesh_.boundary()[patchI];

        if (patch.type() == wallFvPatch::typeName)
        {
            const polyPatch& pPatch = patch.patch();

            const labelList::subList FaceCells = patch.faceCells();

            // Check cells with face on wall
            forAll(patch, patchFaceI)
            {
                label nNeighbours = wallUtils.getPointNeighbours
                (
                    pPatch,
                    patchFaceI,
                    neighbours
                );

                label minFaceI = -1;

                ypatch[patchFaceI] = wallUtils.smallestDist
                (
                    cellCentres[FaceCells[patchFaceI]],
                    pPatch,
                    nNeighbours,
                    neighbours,
                    minFaceI
                );
            }
        }
        else
        {
            ypatch = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    volScalarField::GeometricBoundaryField
    (
        mesh.boundary(),
        scalarField(mesh.nCells(), 0.0),
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(mesh)
{
    doAll();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Correct for mesh geom/topo changes
void Foam::nearWallDist::correct()
{
    if (mesh_.morphing())
    {
        // Update size of GeometricBoundaryField
        forAll(mesh_.boundary(), patchI)
        {
            operator[](patchI).setSize(mesh_.boundary()[patchI].size());
        }
    }

    doAll();
}


// ************************************************************************* //
