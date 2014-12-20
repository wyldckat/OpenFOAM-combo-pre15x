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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "polyMesh.H"
#include "pyramidPointFaceRef.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMesh::correctFaceDirections()
{
    Info<< "void polyMesh::correctFaceDirections() : "
        << "correcting face directions." << endl;

    const cellList& c = allCells();

    const labelList& own = allOwner();

    const pointField& p = allPoints();

    label nChanged = 0;

    forAll (faces_, faceI)
    {
        // Create the owner - face dot product and check face area
        scalar dir =
            (faces_[faceI].centre(p) - c[own[faceI]].centre(p, faces_))
          & faces_[faceI].normal(p);

        if (dir < VSMALL)
        {
            Info<< "bool polyMesh::correctFaceDirections() "
                << "face " << faceI << " points the wrong way. Correcting."
                << " Owner cell: " << own[faceI] << endl;

            nChanged++;

            faces_[faceI] = faces_[faceI].reverseFace();
        }
    }

    if (nChanged > 0)
    {
        Info << "Changed " << nChanged << " faces" << endl;

        // Clear out data that is out of date
        clearGeom();
    }

    Info<< "void polyMesh::correctFaceDirections() : "
        << "finished correcting face directions." << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

