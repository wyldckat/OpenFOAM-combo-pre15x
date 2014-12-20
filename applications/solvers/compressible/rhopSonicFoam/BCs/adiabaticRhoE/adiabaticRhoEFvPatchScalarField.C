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

#include "adiabaticRhoEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adiabaticRhoEFvPatchScalarField::adiabaticRhoEFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


adiabaticRhoEFvPatchScalarField::adiabaticRhoEFvPatchScalarField
(
    const adiabaticRhoEFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


adiabaticRhoEFvPatchScalarField::adiabaticRhoEFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


adiabaticRhoEFvPatchScalarField::adiabaticRhoEFvPatchScalarField
(
    const adiabaticRhoEFvPatchScalarField& ptpsf,
    const scalarField& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void adiabaticRhoEFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void adiabaticRhoEFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void adiabaticRhoEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& rhop =
        lookupPatchField<volScalarField, scalar>("rho");
    const fvPatchField<vector>& rhoUp =
        lookupPatchField<volVectorField, vector>("rhoU");

    valueFraction() = 
        rhop.snGrad()/(rhop.snGrad() - rhop*this->patch().deltaCoeffs());

    refValue() = 0.5*rhop*magSqr(rhoUp/rhop);
    refGrad() =
    (
        refValue() 
      - (
            0.5*rhop.patchInternalField()*
            magSqr(rhoUp.patchInternalField()/rhop.patchInternalField())
        )
    )*patch().deltaCoeffs();

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void adiabaticRhoEFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, adiabaticRhoEFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
