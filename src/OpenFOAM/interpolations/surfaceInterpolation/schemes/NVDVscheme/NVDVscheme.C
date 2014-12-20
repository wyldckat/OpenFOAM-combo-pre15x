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
    Class to create the weighting-factors based on the NVDV
    (Normalised Variable Diagram).
    The particular differencing scheme class is supplied as a template argument,
    the weight function of which is called by the weight function of this class
    for the internal faces as well as faces of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the face-flux, the cell and face
    gradients (from which the normalised variable distribution may be created)
    and the cell centre distance.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Return the interpolation weighting factors

template<class NVDVweight>
tmp<surfaceScalarField> NVDVscheme<NVDVweight>::weights
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

    volTensorField gradc(fvc::grad(phi));

    surfaceVectorField d =
        mesh().Sf()/(mesh().magSf()*mesh().surfaceInterpolation::deltaCoeffs());

    if (!mesh().orthogonal())
    {
        d -=
            mesh().surfaceInterpolation::correctionVectors()
           /mesh().surfaceInterpolation::deltaCoeffs();
    }

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(weights, face)
    {
        weights[face] =
            this->weight
            (
                weights[face],
                faceFlux_[face],
                phi[owner[face]],
                phi[neighbour[face]],
                gradc[owner[face]],
                gradc[neighbour[face]],
                d[face]
            );
    }


    surfaceScalarField::GeometricBoundaryField&
        bWeights = weightingFactors.boundaryField();

    forAll(bWeights, patchI)
    {
        if (bWeights[patchI].coupled())
        {
            scalarField& pWeights = bWeights[patchI];
            const scalarField& pFaceFlux = faceFlux_.boundaryField()[patchI];
            vectorField pphiP =
                phi.boundaryField()[patchI].patchInternalField();
            vectorField pphiN =
                phi.boundaryField()[patchI].patchNeighbourField();
            tensorField pGradcP =
                gradc.boundaryField()[patchI].patchInternalField();
            tensorField pGradcN =
                gradc.boundaryField()[patchI].patchNeighbourField();
            const vectorField& pD = d.boundaryField()[patchI];

            forAll(pWeights, face)
            {
                pWeights[face] =
                    this->weight
                    (
                        pWeights[face],
                        pFaceFlux[face],
                        pphiP[face],
                        pphiN[face],
                        pGradcP[face],
                        pGradcN[face],
                        pD[face]
                    );
            }
        }
    }

    return tWeightingFactors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
