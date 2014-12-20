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
    Reads surface and applies surface regioning to a mesh. Uses boundaryMesh
    to do the hard work.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "boundaryMesh.H"
#include "repatchPolyMesh.H"
#include "faceSet.H"

using namespace Foam;

void patchFaceSet
(
    repatchPolyMesh& mesh,
    const boundaryMesh& bMesh,
    const word& setName,
    const labelList& nearest
)
{
    faceSet faceLabels(mesh, setName);
    Info<< "Read " << faceLabels.size() << " faces to repatch ..." << endl;

    // Add surface patches to end of normal patches. Make unique by prefixing
    // surf_.

    const PtrList<boundaryPatch>& bPatches = bMesh.patches();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> newPatchPtrList(patches.size() + bPatches.size());

    // Copy old ones
    forAll(patches, i)
    {
        const polyPatch& pp = patches[i];

        newPatchPtrList[i] = pp.clone
        (
            mesh.boundaryMesh(),
            i,
            pp.size(),
            pp.start()
        ).ptr();
    }

    // Store old size of patches (since mesh.changePatches changes size)
    label nPatches = patches.size();

    // Add boundary patches.
    forAll(bPatches, i)
    {
        const boundaryPatch& bp = bPatches[i];

        word newName = "surf_" + bp.name();

        Info<< "Added surface patch " << newName << " at position "
            << patches.size() + i << endl;

        newPatchPtrList[nPatches + i] = polyPatch::New
        (
            bp.physicalType(),
            newName,
            0,                  // size
            mesh.nFaces(),      // start
            nPatches + i,
            mesh.boundaryMesh()
        ).ptr();
    }

    // Actually add new list of patches
    mesh.changePatches(newPatchPtrList);


    label nUnchanged = 0;

    for
    (
        faceSet::const_iterator iter = faceLabels.begin();
        iter != faceLabels.end();
        ++iter
    )
    {
        label faceI = iter.key();

        if (faceI >= mesh.nInternalFaces() || faceI < mesh.nFaces())
        {
            label nearestSurfFaceI = nearest[faceI - mesh.nInternalFaces()];

            if (nearestSurfFaceI == -1)
            {
                // Keep face as is
                nUnchanged++;
            }
            else
            {
                label surfPatchI = bMesh.whichPatch(nearestSurfFaceI);

                mesh.changePatchID(faceI, nPatches + surfPatchI);
            }
        }
    }


    if (nUnchanged > 0)
    {
        WarningIn
        (
            "void patchFaceSet"
            "(repatchPolyMesh& mesh, const boundaryMesh& bMesh,"
            "const word& setName, const labelList& nearest)"
        )   << "There are " << nUnchanged << " faces in faceSet "
            << faceLabels.name() << " that did not get repatched since they"
            << " are too far away from the surface" << endl;
    }

    mesh.repatch();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surface file");
    argList::validOptions.insert("faceSet", "faceSet name");

#   include "setRootCase.H"
#   include "createTime.H"

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


    fileName surfName(args.args()[3]);

    Info<< "Reading surface from " << surfName << " ..." << endl;

    bool readSet = args.options().found("faceSet");
    word setName;

    if (readSet)
    {
        setName = args.options()["faceSet"];

        Info<< "Repatching only the faces in faceSet " << setName
            << " according to nearest surface triangle ..." << endl;
    }
    else
    {
        Info<< "Patching all boundary faces according to nearest surface"
            << " triangle ..." << endl;
    }


    Info<< "Before patching:" << nl
        << "    patch\tsize" << endl;

    forAll(mesh.boundaryMesh(), patchI)
    {
        Info<< "    " << mesh.boundaryMesh()[patchI].name() << '\t'
            << mesh.boundaryMesh()[patchI].size() << endl;
    } 
    Info<< endl;


    boundaryMesh bMesh;

    bMesh.readTriSurface(surfName);

    runTime++;

    Info<< "Time now " << runTime.value() << endl;

    if (readSet)
    {
        // Obtain nearest face in bMesh for each boundary face in mesh that
        // is within 'span' of face.
        // Note: should only determine for subset.
        labelList nearest(bMesh.getNearest(mesh, 10.0));

        patchFaceSet(mesh, bMesh, setName, nearest);
    }
    else
    {
        // Obtain nearest face in bMesh for each boundary face in mesh.
        labelList nearest(bMesh.getNearest(mesh, GREAT));

        // Update mesh with new patches
        bMesh.patchify
        (
            nearest,                // nearest bMesh face
            mesh.boundaryMesh(),    // original boundary faces
            mesh
        );
    }


    Info<< "After patching:" << nl
        << "    patch\tsize" << endl;

    forAll(mesh.boundaryMesh(), patchI)
    {
        Info<< "    " << mesh.boundaryMesh()[patchI].name() << '\t'
            << mesh.boundaryMesh()[patchI].size() << endl;
    } 
    Info<< endl;


    // Write resulting mesh

    mesh.write();
    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
