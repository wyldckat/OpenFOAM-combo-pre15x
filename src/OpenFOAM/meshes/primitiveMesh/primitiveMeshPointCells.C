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
#include "cell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcPointCells() const
{
    // Loop through cells and mark up points

    if (debug)
    {
        Info<< "primitiveMesh::calcPointCells() : "
            << "calculating pointCells"
            << endl;
    }

    // It is an error to attempt to recalculate pointCells
    // if the pointer is already set
    if (pcPtr_)
    {
        FatalErrorIn("primitiveMesh::calcPointCells() const")
            << "pointCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellList& cf = cells();

        // Set up temporary storage
        List<DynamicList<label, cellsPerPoint_> > pc(nPoints());

        forAll (cf, cellI)
        {
            const labelList curPoints = cf[cellI].labels(faces());

            forAll (curPoints, pointI)
            {
                pc[curPoints[pointI]].append(cellI);
            }
        }

        // Copy into a plain list
        pcPtr_ = new labelListList(pc.size());
        labelListList& pointCellAddr = *pcPtr_;

        forAll (pc, pointI)
        {
            pointCellAddr[pointI].transfer(pc[pointI].shrink());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::pointCells() const
{
    if (!pcPtr_)
    {
        calcPointCells();
    }

    return *pcPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
