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

#include "inkJetPressureFvPatchScalarField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inkJetPressureFvPatchScalarField::inkJetPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pBar_(0.0),
    pPrime_(0.0),
    frequency_(0.0)
{}


inkJetPressureFvPatchScalarField::inkJetPressureFvPatchScalarField
(
    const inkJetPressureFvPatchScalarField& ijpsf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ijpsf, p, iF, mapper),
    pBar_(ijpsf.pBar_),
    pPrime_(ijpsf.pPrime_),
    frequency_(ijpsf.frequency_)
{}


inkJetPressureFvPatchScalarField::inkJetPressureFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    pBar_(readScalar(dict.lookup("pBar"))),
    pPrime_(readScalar(dict.lookup("pPrime"))),
    frequency_(readScalar(dict.lookup("frequency")))
{
    evaluate();
}


inkJetPressureFvPatchScalarField::inkJetPressureFvPatchScalarField
(
    const inkJetPressureFvPatchScalarField& ijpsf,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(ijpsf, iF),
    pBar_(ijpsf.pBar_),
    pPrime_(ijpsf.pPrime_),
    frequency_(ijpsf.frequency_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void inkJetPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==
    (
        pBar_
       *(1 + pPrime_*::sin(2*physicalConstant::pi*frequency_
       *db().time().value()))
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


// Write
void inkJetPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("pBar") << pBar_ << token::END_STATEMENT << nl;
    os.writeKeyword("pPrime") << pPrime_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, inkJetPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
