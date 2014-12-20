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
    Tet point field reconstructor.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "geometricTetPointFieldReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate point addressing
labelList geometricTetPointFieldReconstructor::procAddressing
(
    const label procNo
) const
{
    // Allocate the addressing
    labelList addr(procMeshes_[procNo].nPoints(), -1);

    const labelList& pointAddr = pointProcAddressing_[procNo];
    const labelList& cellAddr = cellProcAddressing_[procNo];

    label nAddr = 0;

    // Insert point addressing
    forAll (pointAddr, pointI)
    {
        addr[nAddr] = pointAddr[pointI];
        nAddr++;
    }

#   ifdef FACE_DECOMP
    // Insert face addressing.  Only for face decomposition
    const labelList& faceAddr = faceProcAddressing_[procNo];

    const label faceOffset = mesh_.faceOffset();

    forAll (faceAddr, faceI)
    {
        // Remember to decrement the index by one (turning index)
        addr[nAddr] = faceOffset + mag(faceAddr[faceI]) - 1;
        nAddr++;
    }
#   endif

    // Insert cell addressing
    const label cellOffset = mesh_.cellOffset();

    forAll (cellAddr, cellI)
    {
        addr[nAddr] = cellOffset + cellAddr[cellI];
        nAddr++;
    }

    return addr;
}


labelList geometricTetPointFieldReconstructor::procPatchAddressing
(
    const labelList& procToGlobalAddr,
    const label procNo,
    const label patchNo
) const
{
    labelList addr(procMeshes_[procNo].boundary()[patchNo].size(), -1);

    // Algorithm:
    // Go to the global patch, create a lookup list the size of all
    // points in the mesh and then gather the points for the current
    // patch.
    labelList pointLookup(mesh_.nPoints(), -1);

    const labelList& globalPatchPoints =
        mesh_.boundary()[boundaryProcAddressing_[procNo][patchNo]].meshPoints();

    forAll (globalPatchPoints, pointI)
    {
        pointLookup[globalPatchPoints[pointI]] = pointI;
    }

    // Gather the information
    
    const labelList& procPatchPoints =
        procMeshes_[procNo].boundary()[patchNo].meshPoints();

    forAll (procPatchPoints, pointI)
    {
        addr[pointI] =
            pointLookup[procToGlobalAddr[procPatchPoints[pointI]]];
    }

    if (min(addr) < 0)
    {
        FatalErrorIn
        (
            "labelList geometricTetPointFieldReconstructor::"
            "patchProcAddressing\n"
            "(\n"
            "    const labelList& procToGlobalAddr,\n"
            "    const label procNo,\n"
            "    const label patchNo\n"
            ") const"
        )   << "error in addressing"
            << abort(FatalError);
    }

    return addr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
geometricTetPointFieldReconstructor::geometricTetPointFieldReconstructor
(
    tetPolyMesh& mesh,
    const ptrList<tetPolyMesh>& procMeshes,
    const ptrList<labelIOList>& pointProcAddressing,
    const ptrList<labelIOList>& faceProcAddressing,
    const ptrList<labelIOList>& cellProcAddressing,
    const ptrList<labelIOList>& boundaryProcAddressing
)
:
    mesh_(mesh),
    procMeshes_(procMeshes),
    pointProcAddressing_(pointProcAddressing),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing),
    boundaryProcAddressing_(boundaryProcAddressing)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
