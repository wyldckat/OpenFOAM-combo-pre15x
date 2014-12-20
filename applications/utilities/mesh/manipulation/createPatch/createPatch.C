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
    Utility to create patches out of selected boundary faces. Faces come either
    from existing patches or from a faceSet.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "SortableList.H"
#include "ListSearch.H"
#include "repatchPolyMesh.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    //
    // Read control dictionary
    //

    Info<< "Reading createPatchDict\n" << endl;

    IOdictionary patchDict
    (
        IOobject
        (
            "createPatchDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Info<< "Create repatchPolyMesh\n" << endl;

    repatchPolyMesh mesh
    (
        IOobject
        (
            repatchPolyMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    word patchName(patchDict.lookup("name"));

    if (patches.findPatchID(patchName) != -1)
    {
        FatalErrorIn(args.executable())
            << "Patch " << patchName << " already exists in mesh" << endl
            << "Patches are " << patches.names()
            << exit(FatalError);
    }

    word patchType(patchDict.lookup("type"));

    word sourceType(patchDict.lookup("constructFrom"));

    if (sourceType == "patches")
    {
        wordList patchSources(patchDict.lookup("patches"));

        // Add one patch. Keep all the repatched ones (will have 0 faces in
        // them)
        List<polyPatch*> newPatches(patches.size() + 1);

        label newPatchI = 0;
        label meshFaceI = mesh.nInternalFaces();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            Info<< "Copying patch " << pp.name() << " into position "
                << newPatchI << " with faces starting at " << meshFaceI
                << endl;

            newPatches[newPatchI] =
                pp.clone
                (
                    patches,
                    newPatchI,
                    pp.size(),
                    meshFaceI
                ).ptr();

            meshFaceI += pp.size();
            newPatchI++;
        }

        // Add remaining faces as the new patch.
        Info<< "Adding new patch " << patchName << " into position "
            << newPatchI << " with faces starting at " << meshFaceI
            << " and size " << mesh.nFaces() - meshFaceI
            << endl;

        newPatches[newPatchI] =
            polyPatch::New
            (
                patchType,
                patchName,
                mesh.nFaces() - meshFaceI,
                meshFaceI,
                newPatchI,
                patches
            ).ptr();

        // Repatch faces of the patches.
        forAll(patchSources, i)
        {
            label patchI = patches.findPatchID(patchSources[i]);

            if (patchI == -1)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find source patch " << patchSources[i]
                    << endl << "Valid patch names are " << patches.names()
                    << exit(FatalError);
            }

            const polyPatch& pp = patches[patchI];

            Info<< "Moving faces from " << pp.name() << " starting at "
                << pp.start() << " to patch " << newPatchI
                << " starting at " << meshFaceI << endl;

            forAll(pp, i)
            {
                mesh.changePatchID(pp.start() + i, newPatchI);
            }
        }

        // Actually add new list of patches
        mesh.changePatches(newPatches);
    }
    else if (sourceType == "set")
    {
        word setName(patchDict.lookup("set"));

        faceSet faces(mesh, setName);

        Info<< "Read " << faces.size() << " faces from faceSet " << faces.name()
            << endl;


        // Add one patch
        List<polyPatch*> newPatches(patches.size() + 1);

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            Info<< "Copying patch " << pp.name() << " at position "
                << patchI << endl;

            newPatches[patchI] =
                pp.clone
                (
                    patches,
                    patchI,
                    pp.size(),
                    pp.start()
                ).ptr();
        }

        // Add remaining faces as the new patch.
        label patchI = newPatches.size() - 1;

        Info<< "Adding new patch " << patchName << " at position "
            << patchI << " with faces starting at "
            << mesh.nFaces() - faces.size()
            << " and size " << faces.size()
            << endl;

        newPatches[patchI] =
            polyPatch::New
            (
                patchType,
                patchName,
                0,
                mesh.nFaces(),
                patchI,
                patches
            ).ptr();

        // Sort (since faceSet contains faces in arbitrary order)
        labelList faceLabels(faces.toc());

        SortableList<label> patchFaces(faceLabels);

        forAll(patchFaces, i)
        {
            label faceI = patchFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                FatalErrorIn(args.executable())
                    << "Face " << faceI << " specified in set " << faces.name()
                    << " is not an external face of the mesh." << endl
                    << "This application can only repatch existing boundary"
                    << " faces."
                    << exit(FatalError);
            }

            mesh.changePatchID(patchFaces[i], patchI);
        }

        // Actually add new list of patches
        mesh.changePatches(newPatches);
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Fatal source type " << sourceType << endl
            << "Valid source types are 'patches' 'set'"
            << exit(FatalError);
    }

    runTime++;

    mesh.repatch();

    // Write resulting mesh
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
