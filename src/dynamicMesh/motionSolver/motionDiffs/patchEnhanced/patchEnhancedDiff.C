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
    Patch-enhanced motion diffusion.

\*---------------------------------------------------------------------------*/

#include "patchEnhancedDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchEnhancedDiff, 0);
    addToRunTimeSelectionTable(motionDiff, patchEnhancedDiff, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchEnhancedDiff::patchEnhancedDiff
(
    const tetPolyMesh& tetMesh,
    const dictionary& dict
)
:
    motionDiff(tetMesh),
    patchNames_(dict.lookup("distancePatches"))
{
    const polyMesh& m = tetMesh();
    const polyBoundaryMesh& bdry = m.boundaryMesh();

    forAll (patchNames_, i)
    {
        label pID = bdry.findPatchID(patchNames_[i]);

        reduce(pID, maxOp<label>());

        if (pID == -1)
        {
            FatalErrorIn("linearDiff::linearDiff")
                << "Patch " << patchNames_[i] << " specified in dictionary "
                << dict.name() << " not found" << endl
                << "Valid patches are " << bdry.names() << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchEnhancedDiff::~patchEnhancedDiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchEnhancedDiff::enhance(elementScalarField& g) const
{
    const polyMesh& m = tetMesh()();
    const polyBoundaryMesh& bdry = m.boundaryMesh();

    forAll (patchNames_, i)
    {
        label pID = bdry.findPatchID(patchNames_[i]);

        if (pID > -1)
        {
            // Cannot use patch operations: they are made for point fields
            const labelList::subList fc =
                m.boundaryMesh()[pID].faceCells();

            forAll (fc, fcI)
            {
                g[fc[fcI]] *= 2;
            }
        }
    }
}


Foam::tmp<Foam::elementScalarField> Foam::patchEnhancedDiff::gamma() const
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

    // Enhance next to selected patches
    enhance(td());

    return td;
}


// ************************************************************************* //
