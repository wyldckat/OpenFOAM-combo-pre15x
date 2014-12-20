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

Class
    primitiveMesh

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
//#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCellCentresAndVols() const")
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    // Clear the fields for accumulation
    cellCtrs = vector::zero;
    cellVols = 0.0;

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // first estimate the approximate cell centre as the average of face centres

    vectorField cEst(nCells(), vector::zero);
    scalarField nCellFaces(nCells(), 0.0);

    forAll (own, faceI)
    {
        cEst[own[faceI]] += fCtrs[faceI];
        nCellFaces[own[faceI]] += 1;
    }

    forAll (nei, faceI)
    {
        cEst[nei[faceI]] += fCtrs[faceI];
        nCellFaces[nei[faceI]] += 1;
    }

    cEst /= nCellFaces;

    forAll (own, faceI)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =
            max(fAreas[faceI] & (fCtrs[faceI] - cEst[own[faceI]]), VSMALL);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[faceI] + (1.0/4.0)*cEst[own[faceI]];

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs[own[faceI]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVols[own[faceI]] += pyr3Vol;
    }

    forAll (nei, faceI)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =
            max(fAreas[faceI] & (cEst[nei[faceI]] - fCtrs[faceI]), VSMALL);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[faceI] + (1.0/4.0)*cEst[nei[faceI]];

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs[nei[faceI]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVols[nei[faceI]] += pyr3Vol;
    }

    cellCtrs /= cellVols;
    cellVols *= (1.0/3.0);

    //vectorField p(IFstream("points")());
    //cellCtrs = p;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}


const scalarField& primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
