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

Application
    scalePoints

Description
    Scales the mesh points in the polyMesh directory by a factor supplied as an 
    argument.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "polyMesh.H"
#include "argList.H"
#include "IStringStream.H"
#include "boundBox.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("x scaling factor");
    argList::validArgs.append("y scaling factor");
    argList::validArgs.append("z scaling factor");
 
#   include "setRootCase.H"

    scalar xScaleFactor(readScalar(IStringStream(args.args()[3])()));
    scalar yScaleFactor(readScalar(IStringStream(args.args()[4])()));
    scalar zScaleFactor(readScalar(IStringStream(args.args()[5])()));

#   include "createTime.H"

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    boundBox oldGeom(points);

    scalar scaleFactor =
        mag(xScaleFactor) + mag(yScaleFactor) + mag(zScaleFactor);

    if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
    {
        Info<< "Scaling points by " << xScaleFactor << ", "
            << yScaleFactor << ", "
            << zScaleFactor << endl;

        forAll(points, pointI)
        {
            points[pointI].x() *= xScaleFactor;
            points[pointI].y() *= yScaleFactor;
            points[pointI].z() *= zScaleFactor;
        }

        boundBox newGeom(points);
        Info<< "Original bounding box: min = "
            << oldGeom.min() << " max = " << oldGeom.max() << " meters." << endl
            << "Scaled bounding box:   min = "
            << newGeom.min() << " max = " << newGeom.max() << " meters." << endl
            << endl;

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        Info << "Writing points into directory " << points.path() << nl << endl;
        points.write
        (
            runTime.writeFormat(),
            IOstream::currentVersion,
            runTime.writeCompression()
        );
    }
    else
    {
        Info << "Scaling factor set to one.  No scaling performed." << endl;
    }

    return(0);
}


// ************************************************************************* //
