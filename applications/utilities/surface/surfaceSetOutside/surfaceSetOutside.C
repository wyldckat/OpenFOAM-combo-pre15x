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
    Set normal consistent with respect to a user provided 'outside' point.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "orientedSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("Foam surface file");
    argList::validArgs.append("outsidePoint");
    argList::validArgs.append("output file");
    argList args(argc, argv);

    fileName surfFileName(args.args()[1]);
    Info<< "Reading surface from " << surfFileName << endl;

    point outsidePoint(IStringStream(args.args()[2])());
    Info<< "Outside point " << outsidePoint << endl;

    fileName outFileName(args.args()[3]);
    Info<< "Writing surface to " << outFileName << endl;


    // Load surface
    triSurface surf(surfFileName);

    orientedSurface normalSurf(surf, outsidePoint);

    Info<< "Writing new surface to " << outFileName << endl;

    normalSurf.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
