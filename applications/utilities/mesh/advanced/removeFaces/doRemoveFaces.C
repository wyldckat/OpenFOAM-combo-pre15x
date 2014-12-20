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
    Utility to remove faces (combines cells on both sides).

    Takes faceSet of candidates for removal and writes faceSet with faces that
    will actually be removed. (because e.g. would cause two faces between the
    same cells). See removeFaces in dynamicMesh library for constraints.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "removeFaces.H"
#include "ListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return true if all elems are in the set
bool allInSet(const boolList& set, const labelList& elems)
{
    forAll(elems, elemI)
    {
        if (!set[elems[elemI]])
        {
            return false;
        }
    }
    return true;
}


// Set selected elements of boolList.
void set(const labelList& elems, const bool val, boolList& set)
{
    forAll(elems, i)
    {
        set[elems[i]] = val;
    }
}


// Find edges where all faces are in the candidateSet. If so protect all
// neighbours of these faces.
void filterCandidates(const primitiveMesh& mesh, faceSet& candidateSet)
{
    // Convert set to boolList
    boolList candidate(mesh.nFaces());

    forAllConstIter(faceSet, candidateSet, iter)
    {
        candidate[iter.key()] = true;
    }

    // Check for edges whose all faces are in the set.

    const labelListList& edgeFaces = mesh.edgeFaces();

    label nOldStars = -1;

    while (true)
    {
        label nStars = 0;

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = mesh.edgeFaces()[edgeI];

            if (allInSet(candidate, eFaces))
            {
                // Found edge with all faces around it marked for removal
                // (= edge star). Protect all other edges.
                Info<< "Found star around edge:" << edgeI << " verts:"
                    << mesh.edges()[edgeI] << endl;


                // Protect all faces neighbouring this edge-star.
                forAll(eFaces, i)
                {
                    label faceI = eFaces[i];

                    const labelList& fEdges = mesh.faceEdges()[faceI];

                    forAll(fEdges, fEdgeI)
                    {
                        set(mesh.edgeFaces()[fEdges[fEdgeI]], false, candidate);
                    }
                }

                // Reenable collapse of star
                set(eFaces, true, candidate);

                nStars++;
            }
        }

        Info<< "Number of edge-stars:" << nStars << " (was : " << nOldStars
            << " )" << endl;

        if (nStars == nOldStars)
        {
            break;
        }

        nOldStars = nStars;
    }

    // Convert back to faceSet.
    label oldSize = candidateSet.size();
    candidateSet.clear();
    candidateSet.resize(oldSize);

    forAll(candidate, faceI)
    {
        if (candidate[faceI])
        {
            candidateSet.insert(faceI);
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::noParallel();
    Foam::argList::validArgs.append("faceSet");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    word setName(args.args()[3]);

    // Read faces
    faceSet candidateSet(mesh, setName);

    Info<< "Read " << candidateSet.size() << " faces to remove" << nl
        << endl;

    // Remove pointFaces of faces to be removed.
    filterCandidates(mesh, candidateSet);

    {
        faceSet filteredCandidates(mesh, "filteredCandidates", candidateSet);
        Info<< "Writing filtered faceSet to "
            << filteredCandidates.instance()
              /filteredCandidates.local()
              /filteredCandidates.name()
            << endl;

        filteredCandidates.write();
    }

    Info<< "After filtering : " << candidateSet.size() << " faces to remove"
        << nl << endl;


    labelList candidates(candidateSet.toc());

    // Face removal engine. Use default split-face recognition angle (10deg?)
    removeFaces faceRemover(mesh);

    // Get compatible set of faces and connected sets of cells.
    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    faceRemover.compatibleRemoves
    (
        candidates,
        cellRegion,
        cellRegionMaster,
        facesToRemove
    );

    faceSet compatibleRemoves(mesh, "compatibleRemoves", facesToRemove);

    Info<< "Original faces to be removed:" << candidateSet.size() << nl
        << "New faces to be removed:" << compatibleRemoves.size() << nl
        << endl;

    Info<< "Writing new faces to be removed to faceSet "
        << compatibleRemoves.instance()
          /compatibleRemoves.local()
          /compatibleRemoves.name()
        << endl;

    compatibleRemoves.write();


    // Topo changes container
    polyTopoChange meshMod(mesh);

    // Insert mesh refinement into polyTopoChange.
    faceRemover.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Do all changes
    Info<< "Morphing ..." << endl;

    runTime++;

    autoPtr<mapPolyMesh> morphMap = polyTopoChanger::changeMesh(mesh, meshMod);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Write resulting mesh
    Info << "Writing morphed mesh to time " << runTime.value() << endl;

    mesh.write();


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
