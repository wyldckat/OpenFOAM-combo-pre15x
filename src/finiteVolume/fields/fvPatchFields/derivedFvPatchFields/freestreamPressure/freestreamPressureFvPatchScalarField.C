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

#include "freestreamPressureFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freestreamPressureFvPatchScalarField::freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    zeroGradientFvPatchScalarField(p, iF)
{}


freestreamPressureFvPatchScalarField::freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


freestreamPressureFvPatchScalarField::freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict)
{}


freestreamPressureFvPatchScalarField::freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& wbppsf,
    const scalarField& iF
)
:
    zeroGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field

void freestreamPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const freestreamFvPatchVectorField& Up = 
        refCast<const freestreamFvPatchVectorField>
        (
            lookupPatchField<volVectorField, vector>("U")
        );

    const surfaceScalarField& phi = 
        db().lookupObject<surfaceScalarField>("phi");

    fvPatchField<scalar>& phip =
        const_cast<fvPatchField<scalar>&>
        (
            patchField<surfaceScalarField, scalar>(phi)
        );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        phip = patch().Sf() & Up.freestreamValue();
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>("rho");

        phip = rhop*(patch().Sf() & Up.freestreamValue());
    }
    else
    {
        FatalErrorIn("freestreamPressureFvPatchScalarField::updateCoeffs()")
            << "dimensions of phi are not correct"
            << abort(FatalError);
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const freestreamPressureFvPatchScalarField& ptf)
{
    ptf.write(os);

    os.check
    (
        "operator<<(Ostream&, const freestreamPressureFvPatchScalarField&)"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, freestreamPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
