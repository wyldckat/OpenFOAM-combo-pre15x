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

void primitiveMesh::calcCellCells() const
{
    // Loop through faceCells and mark up neighbours

    if (debug)
    {
        Info<< "primitiveMesh::calcCellCells() : calculating cellCells"
            << endl;
    }

    // It is an error to attempt to recalculate cellCells
    // if the pointer is already set
    if (ccPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCellCells() const")
            << "cellCells already calculated"
            << abort(FatalError);
    }
    else
    {
        List<DynamicList<label, facesPerCell_> > cc(nCells());

        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        forAll (nei, faceI)
        {
            cc[own[faceI]].append(nei[faceI]);
            cc[nei[faceI]].append(own[faceI]);
        }

        // Create the storage
        ccPtr_ = new labelListList(cc.size());
        labelListList& cellCellAddr = *ccPtr_;

        // Pack the addressing
        forAll (cc, cellI)
        {
            cellCellAddr[cellI].transfer(cc[cellI].shrink());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::cellCells() const
{
    if (!ccPtr_)
    {
        calcCellCells();
    }

    return *ccPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
