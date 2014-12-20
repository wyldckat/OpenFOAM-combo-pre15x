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

void primitiveMesh::calcPointEdges() const
{
    // Loop through edges and mark up points

    if (debug)
    {
        Info<< "primitiveMesh::calcPointEdges() : "
            << "calculating pointEdges"
            << endl;
    }

    // It is an error to attempt to recalculate pointEdges
    // if the pointer is already set
    if (pePtr_)
    {
        FatalErrorIn("primitiveMesh::calcPointEdges() const")
            << "pointEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        const edgeList& e = edges();

        // Set up temporary storage
        List<DynamicList<label, edgesPerPoint_> > pe(nPoints());

        forAll (e, edgeI)
        {
            pe[e[edgeI].start()].append(edgeI);
            pe[e[edgeI].end()].append(edgeI);
        }

        // Copy into a plain list
        pePtr_ = new labelListList(pe.size());
        labelListList& pointEdgeAddr = *pePtr_;

        forAll (pe, pointI)
        {
            pointEdgeAddr[pointI].transfer(pe[pointI].shrink());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::pointEdges() const
{
    if (!pePtr_)
    {
        calcPointEdges();
    }

    return *pePtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
