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
    Reads .msh file as written by Gmsh.

    Needs surface elements on mesh to be present and aligned with outside faces
    of the mesh. I.e. if the mesh is hexes, the outside faces need to be quads

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "wallPolyPatch.H"
#include "ListSearch.H"
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

    // Get patch labels of meshF.
    labelList newF(meshF.size());
    forAll(newF, fp)
    {
        newF[fp] = meshPointMap[meshF[fp]];
    }

    // Find faces using first point
    const labelList& pFaces = pp.pointFaces()[newF[0]];

    // Go through all these faces and check if there is one which uses all of
    // newF vertices (in any order ;-)
    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        const face& f = pp[faceI];

        // Count uses of vertices of newF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (findIndex(newF, f[fp]) != -1)
            {
                nMatched++;
            }
        }

        if (nMatched == newF.size())
        {
            return faceI;
        }
    }

    FatalErrorIn("findFace") << "Problem : cannot find face " << meshF
        << " in patch" << abort(FatalError);

    return -1;
}



// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append(".msh file");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName mshName(args.args()[3]);

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
                >> prismPoints[3] >> prismPoints[4] >> prismPoints[5]
                >> prismPoints[0] >> prismPoints[1] >> prismPoints[2];

            renumber(mshToFoam, prismPoints);

            cells[cellI++] = cellShape(prism, prismPoints);

            nPrism++;
        }
        else if (elmType == MSHHEX)
        {
            inFile
                >> hexPoints[0] >> hexPoints[1] >> hexPoints[2]
                >> hexPoints[3] >> hexPoints[4] >> hexPoints[5]
                >> hexPoints[6] >> hexPoints[7];

            renumber(mshToFoam, hexPoints);

            cells[cellI++] = cellShape(hex, hexPoints);

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
        DynamicList<face>& pFaces = patchFaces[patchI];

        forAll(pFaces, i)
        {
            const face& f = pFaces[i];

            // Find face in pp using all vertices of f.
            label patchFaceI = findFace(pp, f);

            label meshFaceI = pp.start() + patchFaceI;

            mesh.changePatchID(meshFaceI, patchI);
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
