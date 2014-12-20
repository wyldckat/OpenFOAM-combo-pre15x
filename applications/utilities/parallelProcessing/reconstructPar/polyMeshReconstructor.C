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
    polyMesh reconstructor.

\*---------------------------------------------------------------------------*/

#include "polyMeshReconstructor.H"
#include "Time.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMeshReconstructor::polyMeshReconstructor
(
    polyMesh& mesh,
    const ptrList<polyMesh>& procMeshes,
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMeshReconstructor::reconstructPoints()
{
    IOobject points1Header
    (
        "points",
        procMeshes_[0].time().timeName()/procMeshes_[0].name(),
        procMeshes_[0],
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    // If the points do not exist in the time directory for processor 0 return
    if (!points1Header.headerOk())
    {
        return;
    }

    // Read the field for all the processors
    ptrList<pointIOField> procsPoints(procMeshes_.size());

    forAll (procMeshes_, procI)
    {
        procsPoints.hook
        (
            new pointIOField
            (
                IOobject
                (
                    "points",
                    procMeshes_[procI].time().timeName()
                        /procMeshes_[procI].name(),
                    procMeshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }


    // Create the new points
    vectorField newPoints(mesh_.nPoints());

    forAll (procMeshes_, procI)
    {
        const vectorField& procPoints = procsPoints[procI];

        // Set the cell values in the reconstructed field

        const labelList& pointProcAddressingI = pointProcAddressing_[procI];

        forAll(pointProcAddressingI, pointI)
        {
            newPoints[pointProcAddressingI[pointI]] = procPoints[pointI];
        }
    }

    mesh_.movePoints(newPoints);
    mesh_.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
