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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
tmp<surfaceScalarField> LimitedScheme<Type, Limiter, LimitFunc>::limiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        new surfaceScalarField
        (
            IOobject
            (
                "limiter",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );
    surfaceScalarField& lim = tLimiter();

    tmp<GeometricField<typename Limiter::phiType, fvPatchField, volMesh> >
        tlPhi = LimitFunc<Type>()(phi);

    const GeometricField<typename Limiter::phiType, fvPatchField, volMesh>&
        lPhi = tlPhi();

    GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>
        gradc(fvc::grad(lPhi));

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    surfaceVectorField d =
        mesh.Sf()/(mesh.magSf()*mesh.surfaceInterpolation::deltaCoeffs());

    if (!mesh.orthogonal())
    {
        d -= mesh.surfaceInterpolation::correctionVectors()
           /mesh.surfaceInterpolation::deltaCoeffs();
    }

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    scalarField& pLim = lim.internalField();

    forAll(pLim, face)
    {
        pLim[face] = Limiter::limiter
        (
            CDweights[face],
            this->faceFlux_[face],
            lPhi[owner[face]],
            lPhi[neighbour[face]],
            gradc[owner[face]],
            gradc[neighbour[face]],
            d[face]
        );
    }

    surfaceScalarField::GeometricBoundaryField& bLim =
        lim.boundaryField();

    forAll(bLim, patchI)
    {
        scalarField& pLim = bLim[patchI];

        if (bLim[patchI].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchI];
            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchI];
            Field<typename Limiter::phiType> plPhiP =
                lPhi.boundaryField()[patchI].patchInternalField();
            Field<typename Limiter::phiType> plPhiN =
                lPhi.boundaryField()[patchI].patchNeighbourField();
            Field<typename Limiter::gradPhiType> pGradcP =
                gradc.boundaryField()[patchI].patchInternalField();
            Field<typename Limiter::gradPhiType> pGradcN =
                gradc.boundaryField()[patchI].patchNeighbourField();
            const vectorField& pD = d.boundaryField()[patchI];

            forAll(pLim, face)
            {
                pLim[face] = Limiter::limiter
                (
                    pCDweights[face],
                    pFaceFlux[face],
                    plPhiP[face],
                    plPhiN[face],
                    pGradcP[face],
                    pGradcN[face],
                    pD[face]
                );
            }
        }
        else
        {
            pLim = 1.0;
        }
    }

    return tLimiter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
