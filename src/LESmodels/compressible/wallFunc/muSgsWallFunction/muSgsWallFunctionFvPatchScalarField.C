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

#include "muSgsWallFunctionFvPatchScalarField.H"
#include "LESmodel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESmodels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const muSgsWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const scalarField& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    Istream& is
)
:
    fixedValueFvPatchScalarField(p, iF, is)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const scalarField& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const muSgsWallFunctionFvPatchScalarField& tppsf,
    const scalarField& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void muSgsWallFunctionFvPatchScalarField::evaluate()
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

    const scalarField& muw = lookupPatchField<volScalarField, scalar>("mu");
    const scalarField& rhow = lookupPatchField<volScalarField, scalar>("rho");
    scalarField& muSgsw = *this;

    scalarField magFaceGradU = mag(U.snGrad());

    forAll(muSgsw, facei)
    {
        scalar magUpara = mag(Up[facei]);

        scalar utau = sqrt
        (
            (muSgsw[facei] + muw[facei])
            *magFaceGradU[facei]/rhow[facei]
        );
                
        if(utau > 0)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = kappa*magUpara/utau;
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - utau/(ry[facei]*muw[facei]/rhow[facei])
                    + magUpara/utau 
                    + 1/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));
                        
                scalar df =
                    - 1.0/(ry[facei]*muw[facei]/rhow[facei])
                    - magUpara/sqr(utau)
                    - 1/E*kUu*fkUu/utau;
                    
                scalar utauNew = utau - f/df;
                err = mag((utau - utauNew)/utau);
                utau = utauNew;

            } while (err > 0.01 && ++iter < 10);

            muSgsw[facei] = rhow[facei]*sqr(utau)/magFaceGradU[facei] - muw[facei];
        }
        else
        {
            muSgsw[facei] = 0;
        }
    }

    //This bit just reports the mean y+ value of wall patches
    //Delete if not required
    const scalarField& wallFaceAreas(patch().magSf());
    const scalar totalPatchArea(sum(wallFaceAreas));
    scalar yPlusMean(0.0);
    scalar muSgsMean(0.0);

    forAll(muSgsw, facei)
    {
        scalar utau= sqrt((muSgsw[facei] + muw[facei])*magFaceGradU[facei]/rhow[facei]);

        yPlusMean = yPlusMean +
                    utau/(ry[facei]*muw[facei]/rhow[facei])
                                        *(wallFaceAreas[facei]/totalPatchArea);
        muSgsMean = muSgsMean +
                    muSgsw[facei]*(wallFaceAreas[facei]/totalPatchArea);
    }

    Info << " Patch no " << patch().index() << " = " << patch().type()
         << "; Patch name: " << patch().name()
         << "; y+ Mean = " << yPlusMean
         << "; muSgs(mean)(wall) = " << muSgsMean << endl;

}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, muSgsWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESmodels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
