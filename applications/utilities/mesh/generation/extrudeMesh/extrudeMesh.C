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

Description
    Extrude mesh from existing patch or from patch read from file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "extrudedMesh.H"
#include "dimensionedTypes.H"
#include "linearNormalExtruder.H"
#include "IFstream.H"
#include "faceMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRoots.H"
#   include "createTimeExtruded.H"


    if (args.options().found("sourceRoot") == args.options().found("surface"))
    {
        FatalErrorIn(args.executable())
            << "Need to specify either -sourceRoot/Case/Patch or -surface"
            << " option to specify the source of the patch to extrude"
            << exit(FatalError);
    }


    autoPtr<extrudedMesh> eMeshPtr(NULL);

    if (args.options().found("sourceRoot"))
    {
        fileName rootDirSource(args.options()["sourceRoot"]);
        fileName caseDirSource(args.options()["sourceCase"]);
        fileName patchName(args.options()["sourcePatch"]);

        Info<< "Extruding patch " << patchName
            << " on mesh " << rootDirSource << ' ' << caseDirSource << nl
            << endl;

        Time runTime
        (
            Time::controlDictName,
            rootDirSource,
            caseDirSource
        );

#       include "createPolyMesh.H"

        label patchID = mesh.boundaryMesh().findPatchID(patchName);

        if (patchID == -1)
        {
            FatalErrorIn(args.executable())
                << "Cannot find patch " << patchName
                << " in the source mesh.\n"
                << "Valid patch names are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchID];

        {
            fileName surfName(patchName + ".sMesh");

            Info<< "Writing patch as surfaceMesh to " << surfName << nl << endl;

            faceMesh fMesh(pp.localFaces(), pp.localPoints());

            OFstream os(surfName);
            os << fMesh << nl;
        }

        eMeshPtr.reset
        (
            new extrudedMesh
            (
                IOobject
                (
                    extrudedMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded
                ),
                pp,
                nLayers,                        // number of layers
                linearNormalExtruder(-thickness) // overall thickness (signed!)
            )
        );
    }
    else
    {
        // Read from surface
        fileName surfName(args.options()["surface"]);

        Info<< "Extruding surfaceMesh read from file " << surfName << nl
            << endl;

        IFstream is(surfName);

        faceMesh fMesh(is);

        Info<< "Read patch from file " << surfName << ':' << nl
            << "    points : " << fMesh.points().size() << nl
            << "    faces  : " << fMesh.size() << nl
            << endl;

        eMeshPtr.reset
        (
            new extrudedMesh
            (
                IOobject
                (
                    extrudedMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded
                ),
                fMesh,
                nLayers,                         // number of layers
                linearNormalExtruder(-thickness) // overall thickness (signed!)
            )
        );        
    }

    const extrudedMesh& eMesh = eMeshPtr();

    eMesh.checkMesh();

    if (!eMesh.write())
    {
        FatalErrorIn(args.executable()) << "Failed writing mesh"
            << exit(FatalError);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
