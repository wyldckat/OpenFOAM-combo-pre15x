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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcEdgeVectors() const
{
    if (debug)
    {
        Info<< "primitiveMesh::calcEdgeVectors() : "
            << "calculating edgeVectors"
            << endl;
    }

    // It is an error to attempt to recalculate edgeVectors
    // if the pointer is already set
    if (edgeVectorsPtr_)
    {
        FatalErrorIn("primitiveMesh::calcEdgeVectors() const")
            << "edgeVectors already calculated"
            << abort(FatalError);
    }

    const edgeList& e = edges();

    const pointField& p = points();

    edgeVectorsPtr_ = new vectorField(nEdges());
    vectorField& v = *edgeVectorsPtr_;

    forAll (e, edgeI)
    {
        v[edgeI] = e[edgeI].vec(p);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& primitiveMesh::edgeVectors() const
{
    if (!edgeVectorsPtr_)
    {
        calcEdgeVectors();
    }

    return *edgeVectorsPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
