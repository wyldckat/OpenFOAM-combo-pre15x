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
    A surface analysis tool which sub-sets the triSurface
    to choose only a part of interest. Based on subsetMesh.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Switch.H"
#include "IOdictionary.H"
#include "boundBox.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surfaceSubsetDict");
    argList::validArgs.append("surface file");
    argList::validArgs.append("output file");
    argList args(argc, argv);

    Info<< "Reading dictionary " << args.args()[1] << " ..." << endl;
    IFstream dictFile(args.args()[1]);
    dictionary meshSubsetDict(dictFile);

    Info<< "Reading surface " << args.args()[2] << " ..." << endl;
    triSurface surf1(args.args()[2]);

    Info<< "Original:" << endl
        << "    triangles   :" << surf1.size() << endl
        << "    edges       :" << surf1.nEdges() << endl
        << "    vertices    :" << surf1.nPoints() << endl
        << "    bounding box:" << boundBox(surf1.localPoints())
        << endl << endl;


    labelList markedPoints
    (
        meshSubsetDict.lookup("localPoints")
    );

    labelList markedEdges
    (
        meshSubsetDict.lookup("edges")
    );

    labelList markedFaces
    (
        meshSubsetDict.lookup("faces")
    );

    pointField markedZone
    (
        meshSubsetDict.lookup("zone")
    );

    if ((markedZone.size() != 0) && (markedZone.size() != 2))
    {
        FatalErrorIn(args.executable())
            << "zone specification should be two points, min and max of "
            << "the boundingbox" << endl
            << "zone:" << markedZone
            << exit(FatalError);
    }
            


    Switch addFaceNeighbours
    (
        meshSubsetDict.lookup("addFaceNeighbours")
    );

    // Mark the cells for the subset

    // Faces to subset
    boolList facesToSubset(surf1.size(), false);


    //
    // pick up faces connected to "localPoints"
    //

    if (markedPoints.size() > 0)
    {
        Info << "Found " << markedPoints.size() << " marked point(s)." << endl;

        // pick up cells sharing the point

        forAll (markedPoints, pointI)
        {
            if
            (
                markedPoints[pointI] < 0
             || markedPoints[pointI] >= surf1.nPoints()
            )
            {
                FatalErrorIn(args.executable())
                    << "localPoint label " << markedPoints[pointI]
                    << "out of range."
                    << " The mesh has got "
                    << surf1.nPoints() << " localPoints."
                    << exit(FatalError);
            }

            const labelList& curFaces =
                surf1.pointFaces()[markedPoints[pointI]];

            forAll (curFaces, i)
            {
                facesToSubset[curFaces[i]] =  true;
            }
        }
    }



    //
    // pick up faces connected to "edges"
    //

    if (markedEdges.size() > 0)
    {
        Info << "Found " << markedEdges.size() << " marked edge(s)." << endl;

        // pick up cells sharing the edge

        forAll (markedEdges, edgeI)
        {
            if
            (
                markedEdges[edgeI] < 0
             || markedEdges[edgeI] >= surf1.nEdges()
            )
            {
                FatalErrorIn(args.executable())
                    << "edge label " << markedEdges[edgeI]
                    << "out of range."
                    << " The mesh has got "
                    << surf1.nEdges() << " edges."
                    << exit(FatalError);
            }

            const labelList& curFaces = surf1.edgeFaces()[markedEdges[edgeI]];

            forAll (curFaces, i)
            {
                facesToSubset[curFaces[i]] =  true;
            }
        }
    }


    //
    // pick up faces with centre inside "zone"
    //

    if (markedZone.size() == 2)
    {
        const point& min = markedZone[0];
        const point& max = markedZone[1];

        Info << "Using zone min:" << min << " max:" << max << endl;

        forAll(surf1, faceI)
        {
            const labelledTri& f = surf1[faceI];
            const point centre = f.centre(surf1.points());

            if
            (
                (centre.x() >= min.x())
             && (centre.y() >= min.y())
             && (centre.z() >= min.z())
             && (centre.x() <= max.x())
             && (centre.y() <= max.y())
             && (centre.z() <= max.z())
            )
            {
                facesToSubset[faceI] = true;
            }
        }
    }

    
    //
    // pick up specified "faces"
    //

    // Number of additional faces picked up because of addFaceNeighbours
    label nFaceNeighbours = 0;

    if (markedFaces.size() > 0)
    {
        Info << "Found " << markedFaces.size() << " marked face(s)." << endl;

        // Check and mark faces to pick up
        forAll (markedFaces, faceI)
        {
            if
            (
                markedFaces[faceI] < 0
             || markedFaces[faceI] >= surf1.size()
            )
            {
                FatalErrorIn(args.executable())
                    << "Face label " << markedFaces[faceI] << "out of range."
                    << " The mesh has got "
                    << surf1.size() << " faces."
                    << exit(FatalError);
            }

            // Mark the face
            facesToSubset[markedFaces[faceI]] = true;

            // mark its neighbours if requested
            if (addFaceNeighbours)
            {
                const labelList& curFaces =
                    surf1.faceFaces()[markedFaces[faceI]];

                forAll (curFaces, i)
                {
                    label faceI = curFaces[i];

                    if (!facesToSubset[faceI])
                    {
                        facesToSubset[faceI] =  true;
                        nFaceNeighbours++;
                    }
                }
            }
        }
    }

    if (addFaceNeighbours)
    {
        Info<< "Added " << nFaceNeighbours
            << " faces because of addFaceNeighbours" << endl;
    }

    // Create subsetted surface
    labelList pointMap;
    labelList faceMap;
    triSurface surf2
    (
        surf1.subsetMesh(facesToSubset, pointMap, faceMap)
    );

    Info<< "Subset:" << endl
        << "    triangles   :" << surf2.size() << endl
        << "    edges       :" << surf2.nEdges() << endl
        << "    vertices    :" << surf2.nPoints() << endl
        << "    bounding box:" << boundBox(surf2.localPoints())
        << endl << endl;

//    // Give region numbers the original face number.
//    forAll(surf2, faceI)
//    {
//        surf2[faceI].region() = faceMap[faceI];
//    }

    fileName outFileName(args.args()[3]);

    Info << "Writing surface to " << outFileName << endl;

    surf2.write(outFileName);

    return 0;
}


// ************************************************************************* //
