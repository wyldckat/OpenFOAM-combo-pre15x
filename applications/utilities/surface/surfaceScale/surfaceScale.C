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
    Scale a surface. Like scalePoints but then for surfaces.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "boundBox.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("surface file");
    argList::validArgs.append("x scaling factor");
    argList::validArgs.append("y scaling factor");
    argList::validArgs.append("z scaling factor");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName surfFileName(args.args()[1]);

    Info<< "Reading surf from " << surfFileName << " ..." << endl;

    scalar xScaleFactor(readScalar(IStringStream(args.args()[2])()));
    scalar yScaleFactor(readScalar(IStringStream(args.args()[3])()));
    scalar zScaleFactor(readScalar(IStringStream(args.args()[4])()));

    Info<< "Scaling surface with factors" << nl
        << "    x : " << xScaleFactor << nl
        << "    y : " << yScaleFactor << nl
        << "    z : " << zScaleFactor << nl
        << endl;

    fileName outFileName(args.args()[5]);

    Info<< "Writing surf to " << outFileName << " ..." << endl;


    triSurface surf1(surfFileName);

    pointField newPoints(surf1.points());

    forAll(newPoints, i)
    {
        point& pt = newPoints[i];

        pt.x() *= xScaleFactor;
        pt.y() *= yScaleFactor;
        pt.z() *= zScaleFactor;
    }

    triSurface surf2(surf1, surf1.patches(), newPoints);

    surf2.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
