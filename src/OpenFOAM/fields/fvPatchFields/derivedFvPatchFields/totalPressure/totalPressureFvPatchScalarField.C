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

#include "totalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(p.size(), 0.0)
{}


totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_, mapper)
{}


totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_("p0", dict, p.size())
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
        fvPatchField<scalar>::operator=(p0_);
    }
}


totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& tppsf,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void totalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        lookupPatchField<volVectorField, vector>("U");

    Field<scalar> Upn = Up & patchMesh().nf();

    if (db().foundObject<volScalarField>("psi"))
    {
        const fvPatchField<scalar>& psi =
            lookupPatchField<volScalarField, scalar>("psi");

        //operator==(p0_/(1.0 - 0.5*psi*::min(Upn, 0.0)*mag(Up)));
        operator==(p0_/(1.0 + 0.5*psi*(1.0 - pos(Upn))*magSqr(Up)));
    }
    else
    {
        //operator==(p0_ + 0.5*::min(Upn, 0.0)*mag(Up));
        operator==(p0_ - 0.5*(1.0 - pos(Upn))*magSqr(Up));
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


// Write
void totalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, totalPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
