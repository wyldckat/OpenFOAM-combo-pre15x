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

#include "nuSgsWallFunctionFvPatchScalarField.H"
#include "LESmodel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESmodels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nuSgsWallFunctionFvPatchScalarField::nuSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


nuSgsWallFunctionFvPatchScalarField::nuSgsWallFunctionFvPatchScalarField
(
    const nuSgsWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


nuSgsWallFunctionFvPatchScalarField::nuSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    Istream& is
)
:
    fixedValueFvPatchScalarField(p, iF, is)
{}


nuSgsWallFunctionFvPatchScalarField::nuSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


nuSgsWallFunctionFvPatchScalarField::nuSgsWallFunctionFvPatchScalarField
(
    const nuSgsWallFunctionFvPatchScalarField& tppsf,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nuSgsWallFunctionFvPatchScalarField::evaluate()
{
    const LESmodel& sgsModel
        = db().lookupObject<LESmodel>("turbulenceProperties");

    scalar kappa = dimensionedScalar(sgsModel.lookup("kappa")).value();

    scalar E = dimensionedScalar
    (
        sgsModel.subDict("wallFunctionCoeffs").lookup("E")
    ).value();

    const scalarField& ry = patch().deltaCoeffs();

    const fvPatchVectorField& U = lookupPatchField<volVectorField, vector>("U");
    vectorField Up = U.patchInternalField();

    const scalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");
    scalarField& nuSgsw = *this;


    scalarField magFaceGradU = mag(U.snGrad());

    forAll(nuSgsw, facei)
    {
        scalar magUpara = mag(Up[facei]);

        scalar utau = sqrt((nuSgsw[facei] + nuw[facei])*magFaceGradU[facei]);
                
        if(utau > 0)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = kappa*magUpara/utau;
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - utau/(ry[facei]*nuw[facei])
                    + magUpara/utau 
                    + 1/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));
                        
                scalar df =
                    - 1.0/(ry[facei]*nuw[facei])
                    - magUpara/sqr(utau)
                    - 1/E*kUu*fkUu/utau;
                    
                scalar utauNew = utau - f/df;
                err = mag((utau - utauNew)/utau);
                utau = utauNew;

            } while (utau > 0 && err > 0.01 && ++iter < 10);

            nuSgsw[facei] = sqr(max(utau, 0))/magFaceGradU[facei] - nuw[facei];
        }
        else
        {
            nuSgsw[facei] = 0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nuSgsWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESmodels
} // End namespace Foam

// ************************************************************************* //
