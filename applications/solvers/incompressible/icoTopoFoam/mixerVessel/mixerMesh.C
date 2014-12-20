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

#include "mixerMesh.H"
#include "Time.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixerMesh::calcMovingMasks() const
{
    if (debug)
    {
        Info<< "void mixerMesh::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void mixerMesh::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = allCells();
    const faceList& f = allFaces();

    const labelList& cellAddr =
        cellZones()[cellZones().findZoneID(movingCellsName_)].addressing();

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    const labelList& innerSliderAddr =
        faceZones()[faceZones().findZoneID(innerSliderName_)].addressing();

    forAll (innerSliderAddr, faceI)
    {
        const face& curFace = f[innerSliderAddr[faceI]];

        forAll (curFace, pointI)
        {
            movingPointsMask[curFace[pointI]] = 1;
        }
    }

    const labelList& outerSliderAddr =
        faceZones()[faceZones().findZoneID(outerSliderName_)].addressing();

    forAll (outerSliderAddr, faceI)
    {
        const face& curFace = f[outerSliderAddr[faceI]];

        forAll (curFace, pointI)
        {
            movingPointsMask[curFace[pointI]] = 0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::mixerMesh::mixerMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    fvMesh(io),
    cs_
    (
        "axisCS",
        vector(0, 0, 0), // origin
        vector(0, 0, 1), // axis
        vector(1, 0, 0)  // x-axis
    ),
    rpm_(readScalar(dict.lookup("rpm"))),
    movingCellsName_("movingCells"),
    innerSliderName_("insideSliderZone"),
    outerSliderName_("outsideSliderZone"),
    movingPointsMaskPtr_(NULL)
{
    Info<< "Mixer mesh:" << nl
        << "    origin: " << cs_.origin() << nl
        << "    axis: " << cs_.axis() << nl
        << "    rpm: " << rpm_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixerMesh::~mixerMesh()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::mixerMesh::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMasks();
    }

    return *movingPointsMaskPtr_;
}


void Foam::mixerMesh::slide()
{
     // Rotational speed needs to be converted from rpm
    movePoints
    (
        cs_.toGlobal
        (
            cs_.toLocal(allPoints())
          + vector(0, rpm_*360.0*time().deltaT().value()/60.0, 0)
            *movingPointsMask()
        )
    );

    updateTopology();

    if (morphing())
    {
        if (debug)
        {
            Info << "Mesh is morphing" << endl;
        }

        deleteDemandDrivenData(movingPointsMaskPtr_);
    }

    movePoints
    (
        cs_.toGlobal
        (
            cs_.toLocal(oldAllPoints())
          + vector(0, rpm_*360.0*time().deltaT().value()/60.0, 0)
            *movingPointsMask()
        )
    );
}


// ************************************************************************* //
