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
    for the internal edges as well as edges of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the edge-flux, the cell and edge
    gradients (from which the normalised variable distribution may be created)
    and the cell centre distance.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

\*---------------------------------------------------------------------------*/

#include "areaFields.H"
#include "edgeFields.H"
#include "facGrad.H"
#include "coupledFaPatchFields.H"
#include "upwindEdgeInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline tmp<areaScalarField> limiter(const areaScalarField& phi)
{
    return phi;
}

inline tmp<areaScalarField> limiter(const areaVectorField& phi)
{
    return magSqr(phi);
}

inline tmp<areaScalarField> limiter(const areaTensorField& phi)
{
    return magSqr(phi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Return the interpolation weighting factors

template<class Type, class NVDweight>
tmp<edgeScalarField> faNVDscheme<Type,NVDweight>::weights
(
    const GeometricField<Type, faPatchField, areaMesh>& phi
) const
{
    const faMesh& mesh = this->mesh();

    tmp<edgeScalarField> tWeightingFactors
    (
        new edgeScalarField(mesh.edgeInterpolation::weights())
    );
    edgeScalarField& weightingFactors = tWeightingFactors();

    scalarField& weights = weightingFactors.internalField();

    tmp<areaScalarField> tvf = limiter(phi);
    const areaScalarField& vf = tvf();

    areaVectorField gradc(fac::grad(vf));

//     edgeVectorField d =
//         mesh.Le()
//        /(mesh.magLe()*mesh.edgeInterpolation::deltaCoeffs());

//     if (!mesh.orthogonal())
//     {
//         d -=
//             mesh.edgeInterpolation::correctionVectors()
//            /mesh.edgeInterpolation::deltaCoeffs();
//     }

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const vectorField& n = mesh.faceAreaNormals().internalField();

    forAll(weights, edge)
    {
        vector d = vector::zero;

        if(edgeFlux_[edge] > 0)
        {
            d = mesh.centres()[neighbour[edge]] - mesh.centres()[owner[edge]];
            d -= n[owner[edge]]*(n[owner[edge]]&d);
            d /= mag(d)/mesh.edgeInterpolation::lPN().internalField()[edge];
        }
        else
        {
            d = mesh.centres()[neighbour[edge]] - mesh.centres()[owner[edge]];
            d -= n[neighbour[edge]]*(n[neighbour[edge]]&d);
            d /= mag(d)/mesh.edgeInterpolation::lPN().internalField()[edge];
        }

        weights[edge] =
            this->weight
            (
                weights[edge],
                edgeFlux_[edge],
                vf[owner[edge]],
                vf[neighbour[edge]],
                gradc[owner[edge]],
                gradc[neighbour[edge]],
                d
            );
    }


    GeometricField<scalar, faPatchField, edgeMesh>::GeometricBoundaryField&
        bWeights = weightingFactors.boundaryField();

    forAll(bWeights, patchI)
    {
        if (bWeights[patchI].coupled())
        {
            dynamic_cast<const coupledFaPatchField<Type>&>
                (gradc.boundaryField()[patchI])
                .initPatchNeighbourField();
        }
        else
        {
            scalarField& pWeights = bWeights[patchI];
            unallocLabelList edgeFaces = mesh.boundary()[patchI].edgeFaces();

            vectorField d = 
                mesh.boundary()[patchI].edgeCentres()
              - mesh.boundary()[patchI].edgeFaceCentres();

            forAll(pWeights, edgeI)
            {
                pWeights[edgeI] = 0.0;
//                 pWeights[edgeI] = 0.65;
//                 pWeights[edgeI] =
//                     this->weight
//                     (
//                         pWeights[edgeI],
//                         edgeFlux_.boundaryField()[patchI][edgeI],
//                         vf[edgeFaces[edgeI]],
//                         vf.boundaryField()[patchI][edgeI],
//                         gradc[edgeFaces[edgeI]],
//                         gradc[edgeFaces[edgeI]],
//                         d[edgeI]
//                     );
            }
        }
    }

    forAll(bWeights, patchI)
    {
        if (bWeights[patchI].coupled())
        {
//             scalarField& pWeights = bWeights[patchI];
//             const scalarField& pEdgeFlux = edgeFlux_.boundaryField()[patchI];
//             scalarField pvfP =
//                 vf.boundaryField()[patchI].patchInternalField();
//             scalarField pvfN =
//                 vf.boundaryField()[patchI].patchNeighbourField();
//             vectorField pGradcP =
//                 gradc.boundaryField()[patchI].patchInternalField();
//             vectorField pGradcN =
//                 gradc.boundaryField()[patchI].patchNeighbourField();
//             const vectorField& pD = d.boundaryField()[patchI];

//             forAll(pWeights, edge)
//             {
//                 pWeights[edge] =
//                     this->weight
//                     (
//                         pWeights[edge],
//                         pEdgeFlux[edge],
//                         pvfP[edge],
//                         pvfN[edge],
//                         pGradcP[edge],
//                         pGradcN[edge],
//                         pD[edge]
//                     );
//             }
        }
    }

    return tWeightingFactors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
