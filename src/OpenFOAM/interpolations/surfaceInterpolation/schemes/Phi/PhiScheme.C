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

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "coupledFvPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Return the interpolation weighting factors

template<class PhiWeight>
tmp<surfaceScalarField> PhiScheme<PhiWeight>::weights
(
    const volVectorField& phi
) const
{
    tmp<surfaceScalarField> tWeightingFactors
    (
        new surfaceScalarField(mesh().surfaceInterpolation::weights())
    );
    surfaceScalarField& weightingFactors = tWeightingFactors();

    scalarField& weights = weightingFactors.internalField();

    const surfaceVectorField& Sf = mesh().Sf();
    const surfaceScalarField& magSf = mesh().magSf();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    tmp<surfaceScalarField> tUflux = faceFlux_;

    if (faceFlux_.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const volScalarField& rho = 
            phi.db().objectRegistry::lookupObject<volScalarField>("rho");
        tUflux = faceFlux_/fvc::interpolate(rho);
    }
    else if (faceFlux_.dimensions() != dimVelocity*dimArea)
    {
        FatalErrorIn("PhiScheme<PhiWeight>::weights(const volVectorField& phi)")
            << "dimensions of faceFlux_ are not correct"
            << abort(FatalError);
    }

    const surfaceScalarField& Uflux = tUflux();

    forAll(weights, face)
    {
        weights[face] =
            this->weight
            (
                weights[face],
                Uflux[face],
                phi[owner[face]],
                phi[neighbour[face]],
                Sf[face],
                magSf[face]
            );
    }


    surfaceScalarField::GeometricBoundaryField&
        bWeights = weightingFactors.boundaryField();

    forAll(bWeights, patchI)
    {
        if (bWeights[patchI].coupled())
        {
            scalarField& pWeights = bWeights[patchI];
            const vectorField& pSf = Sf.boundaryField()[patchI];
            const scalarField& pmagSf = magSf.boundaryField()[patchI];
            const scalarField& pFaceFlux = Uflux.boundaryField()[patchI];
            vectorField pphiP =
                phi.boundaryField()[patchI].patchInternalField();
            vectorField pphiN =
                phi.boundaryField()[patchI].patchNeighbourField();

            forAll(pWeights, face)
            {
                pWeights[face] =
                    this->weight
                    (
                        pWeights[face],
                        pFaceFlux[face],
                        pphiP[face],
                        pphiN[face],
                        pSf[face],
                        pmagSf[face]
                    );
            }
        }
    }

    return tWeightingFactors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
