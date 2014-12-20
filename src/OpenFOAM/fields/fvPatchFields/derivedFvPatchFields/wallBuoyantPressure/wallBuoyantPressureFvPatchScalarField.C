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

#include "error.H"

#include "wallBuoyantPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallBuoyantPressureFvPatchScalarField::wallBuoyantPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


wallBuoyantPressureFvPatchScalarField::wallBuoyantPressureFvPatchScalarField
(
    const wallBuoyantPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


wallBuoyantPressureFvPatchScalarField::wallBuoyantPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary&
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


wallBuoyantPressureFvPatchScalarField::wallBuoyantPressureFvPatchScalarField
(
    const wallBuoyantPressureFvPatchScalarField& wbppsf,
    const scalarField& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void wallBuoyantPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& environmentalProperties
        = db().lookupObject<IOdictionary>("environmentalProperties");

    dimensionedVector g(environmentalProperties.lookup("g"));

    const fvPatchField<scalar>& rho =
        lookupPatchField<volScalarField, scalar>("rho");

    gradient() = -rho.snGrad()*(g.value() & patchMesh().Cf());

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, wallBuoyantPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
