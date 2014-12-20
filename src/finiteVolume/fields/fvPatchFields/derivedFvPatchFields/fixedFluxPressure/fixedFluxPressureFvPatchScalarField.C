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

\*---------------------------------------------------------------------------*/

#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& wbppsf,
    const scalarField& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field

void fixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        lookupPatchField<volVectorField, vector>("U");

    const surfaceScalarField& phi = 
        db().lookupObject<surfaceScalarField>("phi");
    fvPatchField<scalar> phip =
        patchField<surfaceScalarField, scalar>(phi);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>("rho");

        phip /= rhop;
    }

    const fvPatchField<scalar>& rAp =
        lookupPatchField<volScalarField, scalar>("1|A(U)");

    gradient() =
    (
        phip - (patch().Sf() & Up)
    )/patch().magSf()/rAp;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const fixedFluxPressureFvPatchScalarField& ptf)
{
    ptf.write(os);

    os.check
    (
        "operator<<(Ostream&, const fixedFluxPressureFvPatchScalarField&)"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedFluxPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
