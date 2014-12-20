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
    Detects and removes sliver faces by merging the opposite edges.

    Does not check whether cells get too few faces so be careful! Takes minimum
    area and sharp angle as parameters. E.g.

        <root> <case> 1E-15 5

    detects faces < 1E-15m^2 area and angles < 5degrees.

    Also might have to be rerun more than once if sliver faces use points
    on other sliver faces.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "morphMesh.H"
#include "mapPolyMesh.H"
#include "physicalConstants.H"
#include "faceCollapser.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Get indices in face f of the two sharp angles. Returns true if two found,
// false otherwise.
bool getSharpAngles
(
    const pointField& points,
    const scalar minCos,
    const label faceI,
    const face& f,

    label& fpA,
    label& fpB
)
{
    fpA = -1;
    fpB = -1;

    forAll(f, fp0)
    {
        label fp1 = (fp0 + 1) % f.size();
        label fpMin1 = (fp0 == 0 ? f.size()-1 : fp0-1);

        // Get the two vectors fp1-fp0 and fpMin1-fp0
        vector e10 = points[f[fp1]] - points[f[fp0]];
        e10 /= mag(e10) + VSMALL;

        vector eMin10 = points[f[fpMin1]] - points[f[fp0]];
        eMin10 /= mag(eMin10) + VSMALL;

        if ((e10 & eMin10) > minCos)
        {
            if (fpA == -1)
            {
                fpA = fp0;
            }
            else if (fpB == -1)
            {
                fpB = fp0;
            }
            else
            {
                // More than 2 sharp corners (so triangle instead of edge)
                // Cannot collapse to edge.
                Warning
                    << "Not collapsing. face:" << faceI << " vertices:" << f
                    << " has sharp corners at "
                    << fpA << points[fpA]
                    << "  "
                    << fpB << points[fpB]
                    << "  "
                    << fp0 << points[fp0]
                    << endl;

                return false;
            }
        }
    }

    if (fpA == -1 || fpB == -1)
    {
        return false;
    }
    else
    {
        return true;

    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("min area");
    argList::validArgs.append("min angle");

#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Create polyMesh for time = " << runTime.value() << nl << endl;
    morphMesh mesh
    (
        IOobject
        (
            morphMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    scalar minArea(readScalar(IStringStream(args.args()[3])()));
    scalar minAngle(readScalar(IStringStream(args.args()[4])()));

    scalar minCos = Foam::cos(minAngle * physicalConstant::pi/180.0);

    Info<< "Collapsing all faces with small area and two vertices with sharp"
        << " angle" << endl
        << "    area :" << minArea << nl
        << "    angle:" << minAngle << nl
        << "    cos  :" << minCos << nl
        << endl;


    DynamicList<label> faceLabels;
    DynamicList<label> startFp;
    DynamicList<label> endFp;


    // Mark off faces that don't need to be visited anymore.
    boolList visited(mesh.nFaces(), false);

    forAll(mesh.faces(), faceI)
    {
        if (!visited[faceI])
        {
            visited[faceI] = true;

            scalar area = mag(mesh.faceAreas()[faceI]);

            if (area < minArea)
            {
                const face& f = mesh.faces()[faceI];

                label fpA, fpB;

                if (getSharpAngles(mesh.points(), minCos, faceI, f, fpA, fpB))
                {
                    Info<< "Marked for collapsing:" << nl
                        << "    face:" << faceI << nl
                        << "    area:" << area << nl
                        << "    vertices:" << f << nl
                        << "    sharp corner at:" << fpA << nl
                        << "    sharp corner at:" << fpB << nl
                        << endl;

                    faceLabels.append(faceI);
                    startFp.append(fpA);
                    endFp.append(fpB);
                }

                // Protect all neighbours from collapse. Bit excessive since
                // only want to make sure that delete vertices on one face are
                // not the kept ones on another face but am lazy here.
                forAll(f, fp)
                {
                    label pointI = f[fp];

                    const labelList& pFaces = mesh.pointFaces()[pointI];

                    forAll(pFaces, i)
                    {
                        visited[pFaces[i]] = true;
                    }
                }
            }
        }
    }

    // Face collapsing engine
    faceCollapser collapser(mesh);

    // Container for all mesh changes
    polyTopoChange meshMod(mesh);

    // Insert mesh refinement into polyTopoChange.
    collapser.setRefinement
    (
        faceLabels.shrink(),
        startFp.shrink(),
        endFp.shrink(),      
        meshMod
    );

    label nRemoved =
        meshMod.removedFaces().size()
      - meshMod.modifiedFaces().size();

    label nChanged = meshMod.modifiedFaces().size();

    Info<< "Changes:" << nl
        << "    faces removed   : " << nRemoved << nl
        << "    faces changed   : " << nChanged << nl
        << endl;

    if (nRemoved != 0 || nChanged != 0)
    {
        // Do all changes
        Info<< "Morphing ..." << endl;

        runTime++;

        mesh.setMorphTimeIndex(runTime.timeIndex());

        mesh.updateTopology(meshMod);

        // Move mesh (since morphing does not do this)
        mesh.movePoints(mesh.morphMap().preMotionPoints());


        // Write resulting mesh
        Info << "Writing collapsed mesh to time " << runTime.value() << endl;

        mesh.write();

        Info<< nl << "You might want to rerun this application to see if"
            << " there are any sliver cells left" << endl;
    }
    else
    {
        Info<< "Mesh not changed ..." << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
