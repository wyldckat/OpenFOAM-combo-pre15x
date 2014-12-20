/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006 Mark Olesen
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

#include "massFluxInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    massFlux_(0),
    phiName_("phi"),
    rhoName_("rho")
{}

massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const massFluxInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    massFlux_(ptf.massFlux_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}

massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    massFlux_(readScalar(dict.lookup("massFlux"))),
    phiName_("phi"),
    rhoName_("rho")
{
    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    if (dict.found("rho"))
    {
        dict.lookup("rho") >> rhoName_;
    }
}

massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const massFluxInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    massFlux_(ptf.massFlux_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}

massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const massFluxInletVelocityFvPatchVectorField& ptf,
    const vectorField& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    massFlux_(ptf.massFlux_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void massFluxInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar avgU = -massFlux_/gSum(patch().magSf());

    vectorField n = patch().nf();

    const surfaceScalarField& phi = db().lookupObject<surfaceScalarField>
    (
        phiName_
    );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(n*avgU);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar> & rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(n*avgU/rhop);
    }
    else
    {
        FatalErrorIn("massFluxInletVelocityFvPatchVectorField::updateCoeffs()")
            << "dimensions of phi are not correct"
            << abort(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void massFluxInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("massFlux") << massFlux_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    massFluxInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
