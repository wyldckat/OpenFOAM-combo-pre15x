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

Description

\*---------------------------------------------------------------------------*/

#include "gammaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gammaContactAngleFvPatchScalarField::gammaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    theta0_(0.0),
    uTheta_(0.0),
    thetaA_(0.0),
    thetaR_(0.0)
{}


gammaContactAngleFvPatchScalarField::gammaContactAngleFvPatchScalarField
(
    const gammaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaA_)
{}


gammaContactAngleFvPatchScalarField::gammaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF),
    theta0_(readScalar(dict.lookup("theta0"))),
    uTheta_(readScalar(dict.lookup("uTheta"))),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR")))
{
    evaluate();
}


gammaContactAngleFvPatchScalarField::gammaContactAngleFvPatchScalarField
(
    const gammaContactAngleFvPatchScalarField& gcpsf,
    const scalarField& iF
)
:
    zeroGradientFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Write
void gammaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("uTheta") << uTheta_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, gammaContactAngleFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
