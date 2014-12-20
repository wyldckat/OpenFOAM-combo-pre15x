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
    Utility to create patches out of selected boundary faces. Faces come either
    from existing patches or from a faceSet.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "SortableList.H"
#include "ListOps.H"
#include "repatchPolyMesh.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}


label getPatch(const polyBoundaryMesh& patches, const word& patchName)
{
    label patchI = patches.findPatchID(patchName);

    if (patchI == -1)
    {
        FatalErrorIn("createPatch")
            << "Cannot find source patch " << patchName
            << endl << "Valid patch names are " << patches.names()
            << exit(FatalError);
    }

    return patchI;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    //
    // Read control dictionary
    //

    Info<< "Create repatchPolyMesh for time = " << runTime.value() << nl
        << endl;
    repatchPolyMesh mesh
    (
        IOobject
        (
            repatchPolyMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );


    Info<< "Reading createPatchDict\n" << endl;

    PtrList<dictionary> patchSources
    (
        IOdictionary
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
        ).lookup("patches")
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // New patches.
    List<polyPatch*> newPatches(patches.size() + patchSources.size());

    // Copy old patches.
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


    // Add new patches constructed according to dictionary.

    // For every boundary face the new patch (or -1 if nothing changed).
    labelList newPatchID(mesh.nInternalFaces(), -1);

    forAll(patchSources, addedI)
    {
        const dictionary& dict = patchSources[addedI];

        word patchName(dict.lookup("name"));

        if (patches.findPatchID(patchName) != -1)
        {
            FatalErrorIn(args.executable())
                << "Patch " << patchName << " already exists in mesh" << endl
                << "Existing patches are " << patches.names()
                << exit(FatalError);
        }

        word patchType(dict.lookup("type"));

        word sourceType(dict.lookup("constructFrom"));

        Info<< "Adding new patch " << patchName << " of type " << patchType
            << " as patch " << patches.size()+addedI << endl;

        // Add an empty patch.
        newPatches[patches.size()+addedI] =
            polyPatch::New
            (
                patchType,
                patchName,
                0,
                mesh.nFaces(),
                patches.size()+addedI,
                patches
            ).ptr();


        if (sourceType == "patches")
        {
            wordList patchSources(dict.lookup("patches"));

            // Repatch faces of the patches.
            forAll(patchSources, sourceI)
            {
                label patchI = getPatch(patches, patchSources[sourceI]);

                const polyPatch& pp = patches[patchI];

                Info<< "Moving faces from patch " << pp.name()
                    << " to patch " << patches.size()+addedI << endl;

                forAll(pp, i)
                {
                    label faceI = pp.start() + i;

                    if (newPatchID[faceI - mesh.nInternalFaces()] != -1)
                    {
                        FatalErrorIn(args.executable())
                            << "Face " << faceI << " on new patch " << patchName
                            << " has already been marked for repatching to"
                            << " patch "
                            << newPatchID[faceI - mesh.nInternalFaces()]
                            << exit(FatalError);
                    }
                    newPatchID[faceI - mesh.nInternalFaces()] =
                        patches.size()+addedI;
                }
            }
        }
        else if (sourceType == "set")
        {
            word setName(dict.lookup("set"));

            faceSet faces(mesh, setName);

            Info<< "Read " << faces.size() << " faces from faceSet "
                << faces.name() << endl;

            // Sort (since faceSet contains faces in arbitrary order)
            labelList faceLabels(faces.toc());

            SortableList<label> patchFaces(faceLabels);

            forAll(patchFaces, i)
            {
                label faceI = patchFaces[i];

                if (mesh.isInternalFace(faceI))
                {
                    FatalErrorIn(args.executable())
                        << "Face " << faceI << " specified in set "
                        << faces.name()
                        << " is not an external face of the mesh." << endl
                        << "This application can only repatch existing boundary"
                        << " faces." << exit(FatalError);
                }

                if (newPatchID[faceI - mesh.nInternalFaces()] != -1)
                {
                    FatalErrorIn(args.executable())
                        << "Face " << faceI << " on new patch " << patchName
                        << " has already been marked for repatching to"
                        << " patch "
                        << newPatchID[faceI - mesh.nInternalFaces()]
                        << exit(FatalError);
                }
                newPatchID[faceI - mesh.nInternalFaces()] =
                    patches.size()+addedI;
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Invalid source type " << sourceType << endl
                << "Valid source types are 'patches' 'set'" << exit(FatalError);
        }
    }

    // Change patch ids
    forAll(newPatchID, i)
    {
        if (newPatchID[i] != -1)
        {
            label faceI = i + mesh.nInternalFaces();

            mesh.changePatchID(faceI, newPatchID[i]);
        }
    }

    // Add new list of patches
    mesh.changePatches(newPatches);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    runTime++;

    // Do all topology changes in one go
    mesh.repatch();

    // Write resulting mesh
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
