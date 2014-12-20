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

void primitiveMesh::calcEdgeCells() const
{
    if (debug)
    {
        Info<< "primitiveMesh::calcEdgeCells() : "
            << "calculating edgeCells"
            << endl;
    }

    // It is an error to attempt to recalculate edgeCells
    // if the pointer is already set
    if (ecPtr_)
    {
        FatalErrorIn("primitiveMesh::calcEdgeCells() const")
            << "edgeCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const labelListList& ce = cellEdges();

        // Set up temporary storage
        List<DynamicList<label, cellsPerEdge_> > ec(nEdges());

        forAll (ce, cellI)
        {
            const labelList& curEdges = ce[cellI];

            forAll (curEdges, edgeI)
            {
                ec[curEdges[edgeI]].append(cellI);
            }
        }

        // Copy into a plain list
        ecPtr_ = new labelListList(ec.size());
        labelListList& edgeCellAddr = *ecPtr_;

        forAll (ec, edgeI)
        {
            edgeCellAddr[edgeI].transfer(ec[edgeI].shrink());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::edgeCells() const
{
    if (!ecPtr_)
    {
        calcEdgeCells();
    }

    return *ecPtr_;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
