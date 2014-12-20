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
    Reads .msh file as written by Gmsh.

    Needs surface elements on mesh to be present and aligned with outside faces
    of the mesh. I.e. if the mesh is hexes, the outside faces need to be quads

    Note: There is something seriously wrong with the ordering written in the
    .msh file. Normal operation is to use the ordering as described
    in the manual. Use the -autoInvert to invert based on the geometry.
    Not very well tested.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "wallPolyPatch.H"
#include "cellModeller.H"
#include "repatchPolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Element type numbers
static label MSHTRI   = 2;
static label MSHQUAD  = 3;
static label MSHTET   = 4;
static label MSHPYR   = 7;
static label MSHPRISM = 6;
static label MSHHEX   = 5;


void renumber
(
    const Map<label>& mshToFoam,
    labelList& labels
)
{
    forAll(labels, labelI)
    {
        labels[labelI] = mshToFoam[labels[labelI]];
    }
}


// Find face in pp which uses all vertices in meshF (in mesh point labels)
label findFace(const primitivePatch& pp, const labelList& meshF)
{
    const Map<label> meshPointMap = pp.meshPointMap();

    // meshF[0] in pp labels.
    if (!meshPointMap.found(meshF[0]))
    {
        Warning<< "Not using gmsh face " << meshF
            << " since zero vertex is not on boundary of polyMesh" << endl;
        return -1;
    }

    // Find faces using first point
    const labelList& pFaces = pp.pointFaces()[meshPointMap[meshF[0]]];

    // Go through all these faces and check if there is one which uses all of
    // meshF vertices (in any order ;-)
    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        const face& f = pp[faceI];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (findIndex(meshF, f[fp]) != -1)
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return faceI;
        }
    }

    Warning << "Could not match gmsh face " << meshF
        << " to any of the outside mesh faces that share the same zeroth point:"
        << IndirectList<face>(pp, pFaces) << endl;

    return -1;
}


// Determine whether cell is inside-out by checking for any wrong-oriented
// face.
bool correctOrientation(const pointField& points, const cellShape& shape)
{
    // Get centre of shape.
    point cc(shape.centre(points));

    // Get outwards pointing faces.
    faceList faces(shape.faces());

    forAll(faces, i)
    {
        const face& f = faces[i];

        vector n(f.normal(points));

        // Check if vector from any point on face to cc points outwards
        if (((points[f[0]] - cc) & n) < 0)
        {
            // Incorrectly oriented
            return false;
        }
    }

    return true;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append(".msh file");
    argList::validOptions.insert("autoInvert", "");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName mshName(args.args()[3]);

    bool autoInvert = args.options().found("autoInvert");

    IFstream inFile(mshName);

    word tag(inFile);

    if (tag != "$NOD")
    {
        FatalErrorIn(args.executable())
            << "Did not find $NOD tag on line " << inFile.lineNumber()
            << exit(FatalError);
    }

    label nVerts;

    inFile >> nVerts;

    Info<< "Read nVerts:" << nVerts << endl << endl;

    pointField points(nVerts);

    Map<label> mshToFoam(nVerts);

    for (label pointI = 0; pointI < nVerts; pointI++)
    {
        label mshLabel;
        scalar xVal, yVal, zVal;

        inFile >> mshLabel >> xVal >> yVal >> zVal;

        point& pt = points[pointI];

        pt.x() = xVal;
        pt.y() = yVal;
        pt.z() = zVal;

        mshToFoam.insert(mshLabel, pointI);
    }
    inFile >> tag;

    if (tag != "$ENDNOD")
    {
        FatalErrorIn(args.executable())
            << "Did not find $ENDNOD tag on line " << inFile.lineNumber()
            << exit(FatalError);
    }

    inFile >> tag;

    if (tag != "$ELM")
    {
        FatalErrorIn(args.executable())
            << "Did not find $ELM tag on line " << inFile.lineNumber()
            << exit(FatalError);
    }

    label nElems;

    inFile >> nElems;

    Info<< "Read nElems:" << nElems << endl << endl;


    // Storage for all cells. Too big. Shrink later
    cellShapeList cells(nElems);

    label cellI = 0;

    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    face triPoints(3);
    face quadPoints(4);
    labelList tetPoints(4);
    labelList pyrPoints(5);
    labelList prismPoints(6);
    labelList hexPoints(8);

    label nTet = 0;
    label nPyr = 0;
    label nPrism = 0;
    label nHex = 0;


    // From gmsh physical region to Foam patch
    Map<label> regionToPatch;

    // Storage for patch faces.
    List<DynamicList<face> > patchFaces(0);


    for (label elemI = 0; elemI < nElems; elemI++)
    {
        label elmNumber, elmType, regPhys, regElem, nNodes;

        inFile >> elmNumber >> elmType >> regPhys >> regElem >> nNodes;

        // regPhys on surface elements is region number.

        labelList verts(nNodes);

        if (elmType == MSHTRI)
        {
            inFile >> triPoints[0] >> triPoints[1] >> triPoints[2];

            renumber(mshToFoam, triPoints);

            Map<label>::iterator regFnd = regionToPatch.find(regPhys);

            label patchI = -1;
            if (regFnd == regionToPatch.end())
            {
                // New region. Allocate patch for it.
                patchFaces.setSize(patchFaces.size() + 1);

                patchI = patchFaces.size()-1;

                Info<< "Mapping region " << regPhys << " to Foam patch "
                    << patchI << endl;
                regionToPatch.insert(regPhys, patchI);
            }
            else
            {
                // Existing patch for region
                patchI = regFnd();
            }

            // Add triangle to correct patchFaces.
            patchFaces[patchI].append(triPoints);
        }
        else if (elmType == MSHQUAD)
        {
            inFile
                >> quadPoints[0] >> quadPoints[1] >> quadPoints[2]
                >> quadPoints[3];

            renumber(mshToFoam, quadPoints);

            Map<label>::iterator regFnd = regionToPatch.find(regPhys);

            label patchI = -1;
            if (regFnd == regionToPatch.end())
            {
                // New region. Allocate patch for it.
                patchFaces.setSize(patchFaces.size() + 1);

                patchI = patchFaces.size()-1;

                Info<< "Mapping region " << regPhys << " to Foam patch "
                    << patchI << endl;
                regionToPatch.insert(regPhys, patchI);
            }
            else
            {
                // Existing patch for region
                patchI = regFnd();
            }

            // Add quad to correct patchFaces.
            patchFaces[patchI].append(quadPoints);
        }
        else if (elmType == MSHTET)
        {
            inFile
                >> tetPoints[0] >> tetPoints[1] >> tetPoints[2]
                >> tetPoints[3];

            renumber(mshToFoam, tetPoints);

            cells[cellI++] = cellShape(tet, tetPoints);

            nTet++;
        }
        else if (elmType == MSHPYR)
        {
            inFile
                >> pyrPoints[0] >> pyrPoints[1] >> pyrPoints[2]
                >> pyrPoints[3] >> pyrPoints[4];

            renumber(mshToFoam, pyrPoints);

            cells[cellI++] = cellShape(pyr, pyrPoints);

            nPyr++;
        }
        else if (elmType == MSHPRISM)
        {
            inFile
                >> prismPoints[0] >> prismPoints[1] >> prismPoints[2]
                >> prismPoints[3] >> prismPoints[4] >> prismPoints[5];

            renumber(mshToFoam, prismPoints);

            cells[cellI] = cellShape(prism, prismPoints);

            const cellShape& cell = cells[cellI];

            if (autoInvert && !correctOrientation(points, cell))
            {
                Info<< "Inverting prism " << cellI << endl;
                // Reorder prism.
                prismPoints[0] = cell[0];
                prismPoints[1] = cell[2];
                prismPoints[2] = cell[1];
                prismPoints[3] = cell[3];
                prismPoints[4] = cell[4];
                prismPoints[5] = cell[5];

                cells[cellI] = cellShape(prism, prismPoints);
            }

            cellI++;

            nPrism++;
        }
        else if (elmType == MSHHEX)
        {
            inFile
                >> hexPoints[0] >> hexPoints[1]
                >> hexPoints[2] >> hexPoints[3]
                >> hexPoints[4] >> hexPoints[5]
                >> hexPoints[6] >> hexPoints[7];

            renumber(mshToFoam, hexPoints);

            cells[cellI] = cellShape(hex, hexPoints);

            const cellShape& cell = cells[cellI];

            if (autoInvert && !correctOrientation(points, cell))
            {
                Info<< "Inverting hex " << cellI << endl;

                // Reorder hex.
                hexPoints[0] = cell[4];
                hexPoints[1] = cell[5];
                hexPoints[2] = cell[6];
                hexPoints[3] = cell[7];
                hexPoints[4] = cell[0];
                hexPoints[5] = cell[1];
                hexPoints[6] = cell[2];
                hexPoints[7] = cell[3];

                cells[cellI] = cellShape(hex, hexPoints);
            }

            cellI++;

            nHex++;
        }
        else
        {
            // Unhandled element. Read the vertices.
            for (label i = 0; i < nNodes; i++)
            {
                label mshLabel;

                inFile >> mshLabel;
            }
        }
    }
    cells.setSize(cellI);

    forAll(patchFaces, patchI)
    {
        patchFaces[patchI].shrink();
    }


    Info<< "Cells:" << endl
        << "    total:" << cells.size() << endl
        << "    hex  :" << nHex << endl
        << "    prism:" << nPrism << endl
        << "    pyr  :" << nPyr << endl
        << "    tet  :" << nTet << endl
        << endl;


    Info<< "Patches:" << nl
        << "Patch\tSize" << endl;

    forAll(patchFaces, patchI)
    {
        Info<< "    " << patchI << '\t' << patchFaces[patchI].size() << endl;
    }
    Info<< endl;


    // Problem is that the orientation of the patchFaces does not have to
    // be consistent with the outwards orientation of the mesh faces. So
    // we have to construct the mesh in two stages:
    // 1. define mesh with all boundary faces in one patch
    // 2. use the read patchFaces to find the corresponding boundary face
    //    and repatch it.



    // Create correct number of patches
    // (but without any faces in it)
    faceListList boundaryFaces(patchFaces.size());

    wordList boundaryPatchNames(boundaryFaces.size());

    forAll(boundaryPatchNames, patchI)
    {
        boundaryPatchNames[patchI] = word("patch") + name(patchI);
    }

    wordList boundaryPatchTypes(boundaryFaces.size(), polyPatch::typeName);
    word defaultFacesType = polyPatch::typeName;
    wordList boundaryPatchPhysicalTypes
    (
        boundaryFaces.size(),
        polyPatch::typeName
    );

    repatchPolyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        points,
        cells,
        boundaryFaces,
        boundaryPatchNames,
        boundaryPatchTypes,
        defaultFacesType,
        boundaryPatchPhysicalTypes
    );


    // Now use the patchFaces to patch up the outside faces of the mesh.

    // Get the patch for all the outside faces (= default patch added as last)
    const polyPatch& pp = mesh.boundaryMesh()[mesh.boundaryMesh().size()-1];

    // Go through all the patchFaces and find corresponding face in pp.
    forAll(patchFaces, patchI)
    {
        const DynamicList<face>& pFaces = patchFaces[patchI];

        forAll(pFaces, i)
        {
            const face& f = pFaces[i];

            // Find face in pp using all vertices of f.
            label patchFaceI = findFace(pp, f);

            if (patchFaceI != -1)
            {
                label meshFaceI = pp.start() + patchFaceI;

                mesh.changePatchID(meshFaceI, patchI);
            }
        }
    }

    //Get polyMesh to write to constant
    runTime.setTime(instant(runTime.constant()), 0);

    mesh.repatch();

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

