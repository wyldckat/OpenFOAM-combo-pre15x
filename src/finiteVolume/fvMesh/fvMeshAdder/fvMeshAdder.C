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

Class
    fvMeshAdder

Description

\*----------------------------------------------------------------------------*/

#include "fvMeshAdder.H"
#include "mapAddedPolyMesh.H"
#include "SortableList.H"
//#include "IOobject.H"
#include "faceCoupleInfo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::fvMeshAdder::fvMeshAdder()
:
    polyMeshAdder()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Inplace add mesh1 to mesh0
void Foam::fvMeshAdder::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo
)
{
    // Save old patch starts (needed in mapping later on)
    labelList patchStarts(mesh0.boundaryMesh().size(), -1);
    forAll(mesh0.boundaryMesh(), patchI)
    {
        patchStarts[patchI] = mesh0.boundaryMesh()[patchI].start();
    }


    // Do all changes to mesh0
    List<polyPatch*> allPatches
    (
        addWithoutPatches
        (
            mesh0,
            mesh1,
            coupleInfo
        )
    );

    // Add patches to new mesh.
    mesh0.removeFvBoundary();
    mesh0.addFvPatches(allPatches);

    // Do the mapping of the stored fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MapVolFields<scalar>(patchStarts, mesh0, mesh1);
    MapVolFields<vector>(patchStarts, mesh0, mesh1);
    MapVolFields<tensor>(patchStarts, mesh0, mesh1);
}


// ************************************************************************* //
