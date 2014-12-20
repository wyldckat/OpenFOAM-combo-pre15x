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

#include "error.H"

#include "primitiveMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcPointFaces() const
{
    // Loop through faces and mark up points

    if (debug)
    {
        Info<< "primitiveMesh::calcPointFaces() : calculating pointFaces"
            << endl;
    }

    // It is an error to attempt to recalculate pointFaces
    // if the pointer is already set
    if (pfPtr_)
    {
        FatalErrorIn("primitiveMesh::calcPointFaces()")
            << "pointFaces already calculated"
            << abort(FatalError);
    }
    else
    {
        const faceList& f = faces();

        // Set up temporary storage
        List<DynamicList<label, facesPerPoint_> > pf(nPoints());

        forAll (f, faceI)
        {
            const labelList& curPoints = f[faceI];

            forAll (curPoints, pointI)
            {
                label ptI = curPoints[pointI];

                if (ptI < 0 || ptI >= nPoints())
                {
                    FatalErrorIn("primitiveMesh::calcPointFaces()")
                        << "Face " << faceI
                        << " contains vertex labels out of range: "
                        << curPoints << " Max point index = " << nPoints()
                        << abort(FatalError);
                }

                pf[ptI].append(faceI);
            }
        }

        // Copy into a plain list
        pfPtr_ = new labelListList(pf.size());
        labelListList& pointFaceAddr = *pfPtr_;

        forAll (pf, pointI)
        {
            pointFaceAddr[pointI].transfer(pf[pointI].shrink());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::pointFaces() const
{
    if (!pfPtr_)
    {
        calcPointFaces();
    }

    return *pfPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
