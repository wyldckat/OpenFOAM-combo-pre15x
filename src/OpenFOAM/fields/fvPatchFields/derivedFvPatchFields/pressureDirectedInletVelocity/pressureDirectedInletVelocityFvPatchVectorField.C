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

#include "pressureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletDir_(p.size())
{}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    inletDir_(ptf.inletDir_, mapper)
{}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& pivpvf,
    const vectorField& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void pressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi = 
        db().lookupObject<surfaceScalarField>("phi");

    const fvPatchField<scalar>& phip =
        patchField<surfaceScalarField, scalar>(phi);

    const vectorField& n = patch().nf();
    scalarField ndmagS = (n & inletDir_)*patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(inletDir_*phip/ndmagS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>("rho");

        operator==(inletDir_*phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorIn
        (
            "pressureDirectedInletVelocityFvPatchVectorField::evaluate()"
        )   << "dimensions of phi are not correct"
            << abort(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void pressureDirectedInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureDirectedInletVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
