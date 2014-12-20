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

#include "inkJetVelocityFvPatchVectorField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inkJetVelocityFvPatchVectorField::inkJetVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    UBar_(vector::zero),
    UPrime_(vector::zero),
    frequency_(0.0)
{}


inkJetVelocityFvPatchVectorField::inkJetVelocityFvPatchVectorField
(
    const inkJetVelocityFvPatchVectorField& ijpsf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ijpsf, p, iF, mapper),
    UBar_(ijpsf.UBar_),
    UPrime_(ijpsf.UPrime_),
    frequency_(ijpsf.frequency_)
{}


inkJetVelocityFvPatchVectorField::inkJetVelocityFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    UBar_(dict.lookup("UBar")),
    UPrime_(dict.lookup("UPrime")),
    frequency_(readScalar(dict.lookup("frequency")))
{
    evaluate();
}


inkJetVelocityFvPatchVectorField::inkJetVelocityFvPatchVectorField
(
    const inkJetVelocityFvPatchVectorField& ijpsf,
    const vectorField& iF
)
:
    fixedValueFvPatchVectorField(ijpsf, iF),
    UBar_(ijpsf.UBar_),
    UPrime_(ijpsf.UPrime_),
    frequency_(ijpsf.frequency_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void inkJetVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==
    (
        scale
        (
            UBar_,
            vector::one
          + UPrime_*::sin(2*physicalConstant::pi*frequency_*db().time().value())
        )
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void inkJetVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("UBar") << UBar_ << token::END_STATEMENT << nl;
    os.writeKeyword("UPrime") << UPrime_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, inkJetVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
