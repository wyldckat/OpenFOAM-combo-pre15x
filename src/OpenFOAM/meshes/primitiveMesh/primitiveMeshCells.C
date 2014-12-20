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

void primitiveMesh::calcCells() const
{
    // Loop through faceCells and mark up neighbours

    if (debug)
    {
        Info<< "primitiveMesh::calcCells() : calculating cells"
            << endl;
    }

    // It is an error to attempt to recalculate cells
    // if the pointer is already set
    if (cfPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCells() const")
            << "cells already calculated"
            << abort(FatalError);
    }
    else
    {
        List<DynamicList<label, facesPerCell_> > cf(nCells());

        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        forAll (own, faceI)
        {
            cf[own[faceI]].append(faceI);
        }

        forAll (nei, faceI)
        {
            cf[nei[faceI]].append(faceI);
        }

        // Create the storage
        cfPtr_ = new cellList(cf.size());
        cellList& cellFaceAddr = *cfPtr_;

        // Pack the addressing
        forAll (cf, cellI)
        {
            cellFaceAddr[cellI].transfer(cf[cellI].shrink());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const cellList& primitiveMesh::cells() const
{
    if (!cfPtr_)
    {
        calcCells();
    }

    return *cfPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
