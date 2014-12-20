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
    Class calculates cell quality measures.

\*---------------------------------------------------------------------------*/

#include "cellQuality.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::cellQuality::cellQuality(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::cellQuality::nonOrthogonality() const
{
    tmp<scalarField> tresult
    (
        new scalarField
        (
            mesh_.nCells(), 0.0
        )
    );
    scalarField& result = tresult();

    scalarField sumArea(mesh_.nCells(), 0.0);

    const vectorField& centres = mesh_.cellCentres();
    const vectorField& areas = mesh_.faceAreas();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll (nei, faceI)
    {
        vector d = centres[nei[faceI]] - centres[own[faceI]];
        vector s = areas[faceI];
        scalar magS = mag(s);

        scalar cosDDotS =
            Foam::acos((d & s)/(mag(d)*magS + VSMALL))
            *180.0/physicalConstant::pi;

        result[own[faceI]] += cosDDotS*magS;
        sumArea[own[faceI]] += magS;

        result[nei[faceI]] += cosDDotS*magS;
        sumArea[nei[faceI]] += magS;
    }

    result /= sumArea;

    return tresult;
}


Foam::tmp<Foam::scalarField> Foam::cellQuality::skewness() const
{
    tmp<scalarField> tresult
    (
        new scalarField
        (
            mesh_.nCells(), 0.0
        )
    );
    scalarField& result = tresult();

    scalarField sumArea(mesh_.nCells(), 0.0);

    const vectorField& cellCtrs = mesh_.cellCentres();
    const vectorField& faceCtrs = mesh_.faceCentres();
    const vectorField& areas = mesh_.faceAreas();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll (nei, faceI)
    {
        scalar dOwn = mag(faceCtrs[faceI] - cellCtrs[own[faceI]]);
        scalar dNei = mag(faceCtrs[faceI] - cellCtrs[nei[faceI]]);

        point faceIntersection =
            cellCtrs[own[faceI]]*dNei/(dOwn+dNei)
          + cellCtrs[nei[faceI]]*dOwn/(dOwn+dNei);

        scalar skewness = 
            mag(faceCtrs[faceI] - faceIntersection)
            /(mag(cellCtrs[nei[faceI]] - cellCtrs[own[faceI]]) + VSMALL);

        scalar magS = mag(areas[faceI]);

        result[own[faceI]] += skewness*magS;
        sumArea[own[faceI]] += magS;

        result[nei[faceI]] += skewness*magS;
        sumArea[nei[faceI]] += magS;
    }

    result /= sumArea;

    return tresult;
}


// ************************************************************************* //
