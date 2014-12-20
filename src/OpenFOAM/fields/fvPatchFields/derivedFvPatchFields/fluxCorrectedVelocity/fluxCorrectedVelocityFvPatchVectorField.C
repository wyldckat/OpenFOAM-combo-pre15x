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

#include "fluxCorrectedVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fluxCorrectedVelocityFvPatchVectorField::fluxCorrectedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{}


fluxCorrectedVelocityFvPatchVectorField::fluxCorrectedVelocityFvPatchVectorField
(
    const fluxCorrectedVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


fluxCorrectedVelocityFvPatchVectorField::fluxCorrectedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{
    fvPatchVectorField::operator=(patchInternalField());

    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    if (dict.found("rho"))
    {
        dict.lookup("rho") >> rhoName_;
    }
}


fluxCorrectedVelocityFvPatchVectorField::fluxCorrectedVelocityFvPatchVectorField
(
    const fluxCorrectedVelocityFvPatchVectorField& fcvpvf,
    const vectorField& iF
)
:
    zeroGradientFvPatchVectorField(fcvpvf, iF),
    phiName_(fcvpvf.phiName_),
    rhoName_(fcvpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void fluxCorrectedVelocityFvPatchVectorField::evaluate()
{
    if (!updated())
    {
        updateCoeffs();
    }

    zeroGradientFvPatchVectorField::evaluate();

    const surfaceScalarField& phi = db().lookupObject<surfaceScalarField>
    (
        phiName_
    );

    const fvPatchField<scalar>& phip =
        patchField<surfaceScalarField, scalar>(phi);

    const vectorField& n = patchMesh().nf();
    const Field<scalar>& magS = patchMesh().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        //operator==(pos(phip)*(*this - n*(n & *this)) + n*phip/magSp);
        operator==(*this - n*(n & *this) + n*phip/magS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(*this - n*(n & *this) + n*phip/(rhop*magS));
    }
    else
    {
        FatalErrorIn("fluxCorrectedVelocityFvPatchVectorField::evaluate()")
            << "dimensions of phi are not correct"
            << abort(FatalError);
    }
}


// Write
void fluxCorrectedVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi")
        << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho")
        << rhoName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, fluxCorrectedVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
