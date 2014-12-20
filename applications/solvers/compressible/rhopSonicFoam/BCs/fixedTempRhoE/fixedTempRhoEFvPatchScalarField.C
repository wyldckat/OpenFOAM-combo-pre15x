/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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

#include "fixedTempRhoEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedTempRhoEFvPatchScalarField::fixedTempRhoEFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Twall_(p.size(), 0.0)
{}


fixedTempRhoEFvPatchScalarField::fixedTempRhoEFvPatchScalarField
(
    const fixedTempRhoEFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Twall_(ptf.Twall_)
{}


fixedTempRhoEFvPatchScalarField::fixedTempRhoEFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    Twall_("Twall", dict, p.size())
{}


fixedTempRhoEFvPatchScalarField::fixedTempRhoEFvPatchScalarField
(
    const fixedTempRhoEFvPatchScalarField& tppsf,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    Twall_(tppsf.Twall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void fixedTempRhoEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& thermodynamicProperties = db().lookupObject<IOdictionary>
    (
        "thermodynamicProperties"
    );

    dimensionedScalar Cv(thermodynamicProperties.lookup("Cv"));

    const fvPatchField<scalar>& rhop =
        lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchVectorField& rhoUp =
        lookupPatchField<volVectorField, vector>("rhoU");

    operator==(rhop*(Cv.value()*Twall_ + 0.5*magSqr(rhoUp/rhop)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


// Write
void fixedTempRhoEFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    Twall_.writeEntry("Twall", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedTempRhoEFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
