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
    Class to create the weighting-factors based on the C-NVD
    (Courant number limited Normalised Variable Diagram).
    The particular differencing scheme class is supplied as a template argument,
    the weight function of which is called by the weight function of this class
    for the internal faces as well as faces of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the face-flux, the cell variable
    values, the cell gradients the face Courant number and the cell centre
    distance.  This version of the CNVD does not use the face gradient in the
    creation of the normalised variables but used the cell-centre variable
    values directly.  This allows a limit to be applied to the estimated
    upwind-cell values as required by the CICSAM scheme of Onno Ubbink.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline tmp<volScalarField> limiter(const volScalarField& phi)
{
    return phi;
}

inline tmp<volScalarField> limiter(const volVectorField& phi)
{
    return magSqr(phi);
}

inline tmp<volScalarField> limiter(const volTensorField& phi)
{
    return magSqr(phi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Return the interpolation weighting factors

template<class Type, class CNVDweight>
tmp<surfaceScalarField> CNVDscheme<Type,CNVDweight>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tWeightingFactors
    (
        new surfaceScalarField(mesh.surfaceInterpolation::weights())
    );
    surfaceScalarField& weightingFactors = tWeightingFactors();

    scalarField& weights = weightingFactors.internalField();

    tmp<volScalarField> tvf = limiter(phi);
    const volScalarField& vf = tvf();

    volVectorField gradc(fvc::grad(vf));

    surfaceScalarField Cof =
        mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux_).interpolate
        (
            fvc::surfaceIntegrate(faceFlux_)
        );

    surfaceVectorField d =
        mesh.Sf()/(mesh.magSf()*mesh.surfaceInterpolation::deltaCoeffs());

    if (!mesh.orthogonal())
    {
        d -= mesh.surfaceInterpolation::correctionVectors()
            /mesh.surfaceInterpolation::deltaCoeffs();
    }

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(weights, face)
    {
        weights[face] =
            this->weight
            (
                weights[face],
                faceFlux_[face],
                vf[owner[face]],
                vf[neighbour[face]],
                gradc[owner[face]],
                gradc[neighbour[face]],
                Cof[face],
                d[face]
            );
    }


    GeometricField<scalar, fvPatchField, surfaceMesh>::GeometricBoundaryField&
        bWeights = weightingFactors.boundaryField();

    forAll(bWeights, patchI)
    {
        if (bWeights[patchI].coupled())
        {
            scalarField& pWeights = bWeights[patchI];
            const scalarField& pFaceFlux = faceFlux_.boundaryField()[patchI];
            scalarField pvfP =
                vf.boundaryField()[patchI].patchInternalField();
            scalarField pvfN =
                vf.boundaryField()[patchI].patchNeighbourField();
            vectorField pGradcP =
                gradc.boundaryField()[patchI].patchInternalField();
            vectorField pGradcN =
                gradc.boundaryField()[patchI].patchNeighbourField();
            const scalarField& pCof = Cof.boundaryField()[patchI];
            const vectorField& pD = d.boundaryField()[patchI];

            forAll(pWeights, face)
            {
                pWeights[face] =
                    this->weight
                    (
                        pWeights[face],
                        pFaceFlux[face],
                        pvfP[face],
                        pvfN[face],
                        pGradcP[face],
                        pGradcN[face],
                        pCof[face],
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
