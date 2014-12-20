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

Application
    createSlider

Description
    Create zones for the mixer vessel case

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionSplit.H"
#include "slidingInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    Info << "Adding zones and mesh modifier to the mixer mesh" << endl;

    // Add an empty zone for cut points
    List<pointZone*> pz(1);

    pz[0] = new pointZone
    (
        "cutPointZone",
        labelList(0),
        0,
        mesh.pointZones()
    );

    // Do face zones for slider
    List<faceZone*> fz(3);

    // Inner slider
    const word innerSliderName("insideSlider");
    const polyPatch& innerSlider =
        mesh.boundaryMesh()[mesh.boundaryMesh().findPatchID(innerSliderName)];

    labelList isf(innerSlider.size());

    forAll (isf, i)
    {
        isf[i] = innerSlider.start() + i;
    }

    fz[0] = new faceZone
    (
        "insideSliderZone",
        isf,
        boolList(innerSlider.size(), false),
        0,
        mesh.faceZones()
    );

    // Outer slider
    const word outerSliderName("outsideSlider");
    const polyPatch& outerSlider =
        mesh.boundaryMesh()[mesh.boundaryMesh().findPatchID(outerSliderName)];

    labelList osf(outerSlider.size());

    forAll (osf, i)
    {
        osf[i] = outerSlider.start() + i;
    }

    fz[1] = new faceZone
    (
        "outsideSliderZone",
        osf,
        boolList(outerSlider.size(), false),
        1,
        mesh.faceZones()
    );

    // Add empty zone for cut faces
    fz[2] = new faceZone
    (
        "cutFaceZone",
        labelList(0),
        boolList(0, false),
        2,
        mesh.faceZones()
    );

    // Cell zone contains moving cells
    List<cellZone*> cz(1);

    regionSplit rs(mesh);

    labelList cellMask(rs.seedMask(mesh.findNearestCell(vector::zero)));

    labelList movingCells(mesh.nCells());
    label nMovingCells = 0;

    forAll (cellMask, cellI)
    {
        if (cellMask[cellI] > 0)
        {
            movingCells[nMovingCells] = cellI;
            nMovingCells++;
        }
    }
    movingCells.setSize(nMovingCells);
    Info << "Number of cells in the moving region: " << nMovingCells << endl;

    cz[0] = new cellZone
    (
        "movingCells",
        movingCells,
        0,
        mesh.cellZones()
    );

    Info << "Adding point and face zones" << endl;
    mesh.addZones(pz, fz, cz);

    // Add a topology modifier

    List<polyMeshModifier*> tm(1);

    tm[0] = new slidingInterface
    (
        "mixerSlider",
        0,
        mesh.morphEngine(),
        "outsideSliderZone",
        "insideSliderZone",
        "cutPointZone",
        "cutFaceZone",
        "outsideSlider",
        "insideSlider",
        slidingInterface::INTEGRAL        // Always integral
    );

    Info << "Adding topology modifiers" << endl;
    mesh.addTopologyModifiers(tm);

    mesh.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
