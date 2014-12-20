// The FOAM Project // File: meshTriangulation.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   meshTriangulation
   \\  /           | Family: mesh
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR
    Mattijs Janssens.

-------------------------------------------------------------------------------
*/

#include "meshTriangulation.H"
#include "polyMesh.H"
#include "faceTriangulation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshTriangulation::isInternalFace
(
    const primitiveMesh& mesh,
    const boolList& includedCell,
    const label faceI
)
{
    if (mesh.isInternalFace(faceI))
    {
        label own = mesh.faceOwner()[faceI];
        label nei = mesh.faceNeighbour()[faceI];

        if (includedCell[own] && includedCell[nei])
        {
            // Neighbouring cell will get included in subset
            // as well so face is internal.
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::meshTriangulation::meshTriangulation()
:
    triSurface(),
    nInternalFaces_(0),
    faceMap_()
{}


// Construct from faces of cells
Foam::meshTriangulation::meshTriangulation
(
    const polyMesh& mesh,
    const label internalFacesPatch,
    const boolList& includedCell
)
:
    triSurface(),
    nInternalFaces_(0),
    faceMap_()
{
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    // All faces to be triangulated.     
    boolList faceIsCut(mesh.nFaces(), false);
    label nFaces = 0;
    label nInternalFaces = 0;
    
    forAll(includedCell, cellI)
    {
        // Include faces of cut cells only.
        if (includedCell[cellI])
        {
            const labelList& cFaces = mesh.cells()[cellI];

            forAll(cFaces, i)
            {
                label faceI = cFaces[i];

                if (!faceIsCut[faceI])
                {
                    // First visit of face.
                    nFaces++;
                    faceIsCut[faceI] = true;

                    // See if would become internal or external face
                    if (isInternalFace(mesh, includedCell, faceI))
                    {
                        nInternalFaces++;
                    }
                }
            }
        }
    }

    Info<< "Subset consists of " << nFaces << " faces out of " << mesh.nFaces()
        << " of which " << nInternalFaces << " are internal" << endl;


    // Find upper limit for number of triangles
    // (can be less if triangulation fails)
    label nTotTri = 0;

    forAll(faceIsCut, faceI)
    {
        if (faceIsCut[faceI])
        {
            nTotTri += faces[faceI].nTriangles(points);
        }
    }

    Info<< "nTotTri : " << nTotTri << endl;


    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Storage for all triangles
    List<labelledTri> triangles(nTotTri);
    faceMap_.setSize(nTotTri);

    label triI = 0;

    // Triangulate internal faces
    forAll(faceIsCut, faceI)
    {
        if (faceIsCut[faceI] && isInternalFace(mesh, includedCell, faceI))
        {
            // Face was internal to the mesh and will be 'internal' to
            // the surface.

            // Triangulate face
            faceTriangulation faceTris(points, faces[faceI]);

            if (faceTris.size() == 0)
            {
                WarningIn("meshTriangulation::meshTriangulation")
                    << "Could not find triangulation for face " << faceI
                    << " vertices " << faces[faceI] << " coords "
                    << IndirectList<point>(points, faces[faceI]) << endl;
            }
            else
            {
                // Copy triangles. Make them internalFacesPatch
                forAll(faceTris, i)
                {
                    const triFace& f = faceTris[i];

                    labelledTri& tri = triangles[triI];

                    tri[0] = f[0];
                    tri[1] = f[1];
                    tri[2] = f[2];

                    tri.region() = internalFacesPatch;

                    faceMap_[triI] = faceI;

                    triI++;
                }
            }
        }
    }
    nInternalFaces_ = triI;


    // Triangulate external faces
    forAll(faceIsCut, faceI)
    {
        if (faceIsCut[faceI] && !isInternalFace(mesh, includedCell, faceI))
        {
            // Face will become outside of the surface.

            // Triangulate face
            faceTriangulation faceTris(points, faces[faceI]);

            if (faceTris.size() == 0)
            {
                WarningIn("meshTriangulation::meshTriangulation")
                    << "Could not find triangulation for face " << faceI
                    << " vertices " << faces[faceI] << " coords "
                    << IndirectList<point>(points, faces[faceI]) << endl;
            }
            else
            {
                label patchI = -1;
                bool reverse = false;

                if (mesh.isInternalFace(faceI))
                {
                    patchI = internalFacesPatch;

                    // Check orientation. Check which side of the face gets
                    // included (note: only one side is).
                    if (includedCell[mesh.faceOwner()[faceI]])
                    {
                        reverse = false;
                    }
                    else
                    {
                        reverse = true;
                    }
                }
                else
                {
                    // Face was already outside so orientation ok.

                    patchI = patches.whichPatch(faceI);

                    reverse = false;
                }

                // Copy triangles. Optionally reverse them
                forAll(faceTris, i)
                {
                    const triFace& f = faceTris[i];

                    labelledTri& tri = triangles[triI];

                    if (reverse)
                    {
                        tri[0] = f[0];
                        tri[2] = f[1];
                        tri[1] = f[2];
                    }
                    else
                    {
                        tri[0] = f[0];
                        tri[1] = f[1];
                        tri[2] = f[2];
                    }

                    tri.region() = patchI;

                    faceMap_[triI] = faceI;

                    triI++;
                }
            }
        }
    }

    // Shrink if nessecary (because of invalid triangulations)
    triangles.setSize(triI);
    faceMap_.setSize(triI);

    Info<< "nInternalFaces_:" << nInternalFaces_ << endl;
    Info<< "triangles:" << triangles.size() << endl;


    geometricSurfacePatchList surfPatches(patches.size());

    forAll(patches, patchI)
    {
        surfPatches[patchI] =
            geometricSurfacePatch
            (
                patches[patchI].physicalType(),
                patches[patchI].name(),
                patchI
            );
    }

    // Create globally numbered tri surface
    triSurface globalSurf(triangles, surfPatches, points);

    // Create locally numbered tri surface
    triSurface::operator=
    (
        triSurface
        (
            globalSurf.localFaces(),
            surfPatches,
            globalSurf.localPoints()
        )
    );
}


// ************************************************************************* //
