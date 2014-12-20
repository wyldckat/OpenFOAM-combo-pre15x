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

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "argList.H"
#include "OFstream.H"
#include "surfaceIntersection.H"

using namespace Foam;
using namespace Foam::triSurfaceTools;

labelList countBins
(
    const scalar min,
    const scalar max,
    const label nBins,
    const scalarField& vals
)
{
    scalar dist = nBins/(max - min);

    labelList binCount(nBins, 0);

    forAll(vals, i)
    {
        scalar val = vals[i];

        label index = -1;

        if (Foam::mag(val - min) < SMALL)
        {
            index = 0;
        }
        else if (Foam::mag(val - max) < SMALL)
        {
            index = nBins - 1;
        }
        else
        {
            index = label((vals[i] - min)*dist);

            if ((index < 0) || (index >= nBins))
            {
                Warning
                    << "countBins(const scalar, const scalar, const label"
                    << ", const scalarField&) : "
                    << "value " << vals[i] << " at index " << i
                    << " outside range " << min << " .. " << max << endl;

                if (index < 0)
                {
                    index = 0;
                }
                else
                {
                    index = nBins - 1;
                }
            }
        }
        binCount[index]++;
    }

    return binCount;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    argList::validArgs.clear();
    argList::validOptions.insert("noCleanup", "");
    argList::validArgs.append("surface file");
    argList args(argc, argv);

    fileName surfFileName(args.args()[1]);
    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    //
    // Read and clean
    //

    triSurface surf(surfFileName);

    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    if (args.options().found("noCleanup"))
    {
        Info << "Not cleaning up surface" << endl;
    }
    else
    {
        Info<< "Cleaning surface ..." << endl;
        surf.cleanup(true);
        Info<< endl;

        Info<< "Statistics:" << endl;
        surf.writeStats(Info);
        Info<< endl;
    }

    //
    // Triangle quality
    //

    scalarField triQ(surf.size());
    forAll(surf, faceI)
    {
        const labelledTri& f = surf[faceI];

        triQ[faceI] = triPointRef
        (
            surf.points()[f[0]],
            surf.points()[f[1]],
            surf.points()[f[2]]
        ).quality();
    }
    labelList binCount = countBins(0, 1, 20, triQ);

    Info<< "Triangle quality (equilateral=1, collapsed=0):"
        << endl;


    OSstream& os = Info;
    os.width(4);

    scalar dist = (1.0 - 0.0)/20.0;
    scalar min = 0;
    forAll(binCount, binI)
    {
        Info<< "    " << min << " .. " << min+dist << "  : "
            << 1.0/surf.size() * binCount[binI]
            << endl;
        min += dist; 
    }
    Info<< endl;


    //
    // Check manifold
    //
    DynamicList<label> problemFaces(surf.size()/100 + 1);

    const labelListList& eFaces = surf.edgeFaces();

    label nSingleEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 1)
        {
            problemFaces.append(myFaces[0]);

            nSingleEdges++;
        }
    }
            
    label nMultEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() > 2)
        {
            forAll(myFaces, myFaceI)
            {
                problemFaces.append(myFaces[myFaceI]);
            }

            nMultEdges++;
        }
    }
    problemFaces.shrink();

    if ((nSingleEdges != 0) || (nMultEdges != 0))
    {
        Info<< "Surface is not closed since not all edges connected to "
            << "two faces:" << endl
            << "    connected to one face : " << nSingleEdges << endl
            << "    connected to >2 faces : " << nMultEdges << endl;

        Info<< "Conflicting face labels:" << problemFaces.size() << endl;

        fileName fName("problemFaces");

        Info<< "Dumping conflicting face labels to " << fName << endl
            << "Paste this into the input for surfaceSubset" << endl;

        OFstream fStream(fName);

        fStream << problemFaces;
    }
    else
    {
        Info<< "Surface is closed. All edges connected to two faces." << endl;
    }
    Info<< endl;



    //
    // Check singly connected domain
    //

    labelList faceZone;
    label numZones = surf.markZones(boolList(surf.nEdges(), false), faceZone);

    Info<< "Number of unconnected parts : " << numZones << endl;

    if (numZones > 1)
    {
        Info<< "Splitting surface into parts ..." << endl << endl;

        fileName surfFileNameBase(surfFileName.name());

        for(label zone = 0; zone < numZones; zone++)
        {
            boolList includeMap(surf.size(), false);

            forAll(faceZone, faceI)
            {
                if (faceZone[faceI] == zone)
                {
                    includeMap[faceI] = true;
                }
            }

            labelList pointMap;
            labelList faceMap;

            triSurface subSurf
            (
                surf.subsetMesh
                (
                    includeMap,
                    pointMap,
                    faceMap
                )
            );

            fileName subFileName
            (
                surfFileNameBase.lessExt()
              + "_"
              + name(zone)
              + ".ftr"
            );

            Info<< "writing part " << zone << " size " << subSurf.size()
                << " to " << subFileName << endl;

            subSurf.write(subFileName);
        }

        return 0;
    }



    //
    // Check self-intersection
    //

    triSurfaceSearch querySurf(surf);
    surfaceIntersection inter(querySurf);

    if ((inter.cutEdges().size() == 0) && (inter.cutPoints().size() == 0))
    {
        Info<< "Surface is not self-intersecting" << endl;
    }
    else
    {
        Info<< "Surface is self-intersecting" << endl;
        Info<< "Writing edges of intersection to selfInter.obj" << endl;

        OFstream intStream("selfInter.obj");
        forAll(inter.cutPoints(), cutPointI)
        {
            const point& pt = inter.cutPoints()[cutPointI];

            intStream << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
                << endl;
        }
        forAll(inter.cutEdges(), cutEdgeI)
        {
            const edge& e = inter.cutEdges()[cutEdgeI];

            intStream << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
        }
    }
    Info<< endl;


    // Check orientation
    boolList borderEdge(surf.checkOrientation(false));

    //
    // Colour all faces into zones using borderEdge
    //
    labelList normalZone;
    label numNormalZones = surf.markZones(borderEdge, normalZone);

    Info<< endl
        << "Number of zones (connected area with consistent normal) : "
        << numNormalZones << endl;

    if (numNormalZones > 1)
    {
        Info<< "More than one normal orientation." << endl;
    }
    Info<< endl;
     

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
