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
    Skewness-based motion diffusion.

\*---------------------------------------------------------------------------*/

#include "skewnessDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "elementFields.H"
#include "cellQuality.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(skewnessDiff, 0);
    addToRunTimeSelectionTable(motionDiff, skewnessDiff, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::skewnessDiff::skewnessDiff
(
    const tetPolyMesh& tetMesh,
    const dictionary&
)
:
    motionDiff(tetMesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::skewnessDiff::~skewnessDiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::elementScalarField> Foam::skewnessDiff::gamma() const
{
    tmp<elementScalarField> td
    (
        new elementScalarField 
        (
            IOobject
            (
                "motionDiff",
                tetMesh().time().timeName(),
                tetMesh()(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("1.0", dimless, 1.0)
        )
    );

    scalarField skewness = cellQuality(tetMesh()()).skewness();

    scalar maxSkewness = max(skewness);
    scalar minSkewness = min(skewness);

    // Fiddle factor: ZT?
    scalar p = 0.1;

    scalar delta = p*(maxSkewness - minSkewness);

    if ((maxSkewness - minSkewness) > SMALL)
    {
        td().internalField() +=
            (skewness - minSkewness)
            /(maxSkewness - skewness + delta);
    }

    return td;
}


// ************************************************************************* //
