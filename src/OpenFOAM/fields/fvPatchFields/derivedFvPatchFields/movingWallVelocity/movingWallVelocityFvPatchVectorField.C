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

\*---------------------------------------------------------------------------*/

#include "movingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& pivpvf,
    const vectorField& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void movingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& patch = patchMesh();
    const polyPatch& pp = patch.patch();
    const fvMesh& mesh = patch.boundaryMesh().mesh();


    const pointField& oldAllPoints = mesh.oldAllPoints();

    vectorField oldFc(pp.size());

    forAll(oldFc, i)
    {
        oldFc[i] = pp[i].centre(oldAllPoints);
    }

    vectorField U = (pp.faceCentres() - oldFc)/mesh.time().deltaT().value();


    const scalarField& phi = patchField<surfaceScalarField, scalar>(mesh.phi());

    const vectorField& n = patch.nf();
    const scalarField& magSf = patch.magSf();
    scalarField Un = phi/(magSf + VSMALL);


    vectorField::operator=(U + n*(Un - (n & U)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void movingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, movingWallVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
