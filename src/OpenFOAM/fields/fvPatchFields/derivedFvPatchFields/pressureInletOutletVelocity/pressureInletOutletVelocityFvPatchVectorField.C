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

#include "pressureInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{
    if (ptf.tangentialVelocity_.size())
    {
        tangentialVelocity_ = mapper(ptf.tangentialVelocity_);
    }
}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;

    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    if (dict.found("rho"))
    {
        dict.lookup("rho") >> rhoName_;
    }

    if (dict.found("tangentialVelocity"))
    {
        tangentialVelocity_ = vectorField("tangentialVelocity", dict, p.size());
    }
}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& pivpvf,
    const vectorField& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    tangentialVelocity_(pivpvf.tangentialVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void pressureInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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
        refValue() = n*phip/magS;
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        refValue() = n*phip/(rhop*magS);
    }
    else
    {
        FatalErrorIn
        (
            "pressureInletOutletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << abort(FatalError);
    }

    if (tangentialVelocity_.size())
    {
        refValue() += tangentialVelocity_ - n*(n & tangentialVelocity_);
    }

    valueFraction() = 1.0 - pos(phip);

    mixedFvPatchVectorField::updateCoeffs();
}


// Write
void pressureInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi")
        << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho")
        << rhoName_ << token::END_STATEMENT << nl;
    if (tangentialVelocity_.size())
    {
        tangentialVelocity_.writeEntry("tangentialVelocity", os);
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, pressureInletOutletVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
