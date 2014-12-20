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
    Rotates the mesh so that a patch aligns with normal vector.
    Patch and normal vector are specified in alignPatchDict

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "readAlignPatchDict.H"

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    runTime++;

    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        FatalErrorIn
        (
            "alignPatch.C"
        )   << patchName << "  patch undefined in mesh"
            << exit(FatalError);
    }
    
    vector patchSf = mesh.Sf().boundaryField()[patchID][0];
    patchSf /= mag(patchSf);
    const pointField& oldPoints = mesh.allPoints();
    tensor T = transformationTensor(patchSf, planeNormal);

    pointField newPoints(oldPoints.size());
    newPoints = transform(T, oldPoints);
    mesh.movePoints(newPoints);
    Info << "Writing mesh to time " << runTime.timeName() <<"\n" <<endl;
    mesh.write();

    if (Uheader.headerOk())
    {
        Info<< "Reading U" << endl;
        volVectorField U(Uheader, mesh);
        U == (T & U);

        Info << "Writing U to time " << runTime.timeName() <<"\n" <<endl;
        U.write();
    }

    Info << "End" << endl;

    return 0;
}


// ************************************************************************* //
