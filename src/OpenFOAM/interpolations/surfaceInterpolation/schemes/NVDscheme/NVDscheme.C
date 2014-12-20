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
    Class to create the weighting-factors based on the NVD
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
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
inline tmp<volScalarField> magSqrLimitFunc<Type>::operator()
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    return magSqr(phi);
}


template<>
inline tmp<volScalarField> magSqrLimitFunc<scalar>::operator()
(
    const volScalarField& phi
) const
{
    return phi;
}


template<>
inline tmp<volScalarField> magSqrLimitFunc<tensor>::operator()
(
    const volTensorField& phi
) const
{
    return tr(phi);
}


template<class Type>
inline tmp<volScalarField> rhoMagSqrLimitFunc<Type>::operator()
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::lookupObject<volScalarField>("rho");
    return magSqr(phi/rho);
}


template<>
inline tmp<volScalarField> rhoMagSqrLimitFunc<scalar>::operator()
(
    const volScalarField& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::lookupObject<volScalarField>("rho");
    return phi/rho;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Return the interpolation weighting factors

template<class Type, class NVDweight, template<class> class LimitFunc>
tmp<surfaceScalarField> NVDscheme<Type,NVDweight,LimitFunc>::weights
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

    surfaceVectorField d =
        mesh.Sf()/(mesh.magSf()*mesh.surfaceInterpolation::deltaCoeffs());

    if (!mesh.orthogonal())
    {
        d -=
            mesh.surfaceInterpolation::correctionVectors()
           /mesh.surfaceInterpolation::deltaCoeffs();
    }

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    tmp<volScalarField> tvf = LimitFunc<Type>()(phi);
    const volScalarField& vf = tvf();

    volVectorField gradc(fvc::grad(vf));

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
                        pD[face]
                    );
            }
        }
    }

    /*
    const surfaceScalarField& linearWeights = 
        mesh.surfaceInterpolation::weights();

    surfaceScalarField upwindWeights = upwind<scalar>(vf, faceFlux_).weights();

    surfaceScalarField limiter = 
        (upwindWeights - weightingFactors)/(upwindWeights - linearWeights);

    volScalarField cellLimiter
    (
        IOobject
        (
            "cellLimiter",
            vf.instance(),
            mesh
        ),
        mesh,
        dimless,
        scalarField(vf.size(), 1.0),
        limiter.boundaryField()
    );

    forAll(limiter, i)
    {
        //cellLimiter[owner[i]] = min(cellLimiter[owner[i]], limiter[i]);
        //cellLimiter[neighbour[i]] = min(cellLimiter[neighbour[i]], limiter[i]);

        cellLimiter[owner[i]] *= limiter[i];
        cellLimiter[neighbour[i]] *= limiter[i];
    }

    forAll(limiter, i)
    {
        limiter[i] = min(cellLimiter[owner[i]], cellLimiter[neighbour[i]]);
    }

    weightingFactors = limiter*linearWeights + (1.0 - limiter)*upwindWeights;
    */

    return tWeightingFactors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
