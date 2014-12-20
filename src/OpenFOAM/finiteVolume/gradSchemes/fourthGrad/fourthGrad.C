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

#include "fourthGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
fourthGrad<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf
)
{
    // Algorithm:
    // The fourth-order gradient is calculated in two passes.  First,
    // the standard least-square gradient is assembled.  Then, the
    // fourth-order correction is added to the second-order accurate
    // gradient to complete the accuracy.  However, it is necessary (I
    // think) to temporarily store the second-order gradient so that
    // the higher-order correction is added consistently; the
    // alternative would mean that for some cells the fourth-order
    // correction is added using the mixture of the second-order
    // gradient, incomplete fourth-order gradient etc, which does not
    // sound as a good idea.  

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    // Assemble the second-order least-square gradient

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tfGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "grad("+vsf.name()+')',
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& fGrad = tfGrad();

    // Get reference to least square vectors
    const surfaceVectorField& ownerLs =
        mesh.surfaceInterpolation::leastSquarePvectors();

    const surfaceVectorField& neighbourLs =
        mesh.surfaceInterpolation::leastSquareNvectors();

    // owner/neighbour addressing
    const unallocLabelList& Owner = mesh.owner();
    const unallocLabelList& Neighbour = mesh.neighbour();

    // Internal faces
    forAll(Owner, faceI)
    {
        fGrad[Owner[faceI]] +=
            ownerLs[faceI]*(vsf[Neighbour[faceI]] - vsf[Owner[faceI]]);

        fGrad[Neighbour[faceI]] -=
            neighbourLs[faceI]*(vsf[Neighbour[faceI]] - vsf[Owner[faceI]]);
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchI)
    {
        const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchI];
        const fvPatchVectorField& patchOwnerLs =
            ownerLs.boundaryField()[patchI];
        const labelList::subList FaceCells =
            fGrad.boundaryField()[patchI].patchMesh().faceCells();

        forAll(patchVsf, patchFaceI)
        {
            fGrad[FaceCells[patchFaceI]] +=
                patchOwnerLs[patchFaceI]
               *(patchVsf[patchFaceI] - vsf[FaceCells[patchFaceI]]);
        }
    }

    fGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, fGrad);

    // Assemble the fourth-order gradient

    // Make a copy of the second-order gradient
    GeometricField<GradType, fvPatchField, volMesh> secondfGrad = fGrad;

    // Reconstruct the d-vectors
    surfaceVectorField d =
        mesh.Sf()/(mesh.magSf()*mesh.surfaceInterpolation::deltaCoeffs());

    if (!mesh.orthogonal())
    {
        d -=
            mesh.surfaceInterpolation::correctionVectors()
           /mesh.surfaceInterpolation::deltaCoeffs();
    }


    // Internal faces
    forAll(Owner, faceI)
    {
        fGrad[Owner[faceI]] +=
            0.5*ownerLs[faceI]
           *(
                d[faceI]
              & (secondfGrad[Neighbour[faceI]] - secondfGrad[Owner[faceI]])
            );

        fGrad[Neighbour[faceI]] -=
            0.5*neighbourLs[faceI]
           *(
                d[faceI]
              & (secondfGrad[Neighbour[faceI]] - secondfGrad[Owner[faceI]])
            );
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchI)
    {
        const fvPatchVectorField& patchOwnerLs =
            ownerLs.boundaryField()[patchI];

        const fvPatchVectorField& patchD = d.boundaryField()[patchI];

        const labelList::subList FaceCells =
            fGrad.boundaryField()[patchI].patchMesh().faceCells();

        if (secondfGrad.boundaryField()[patchI].coupled())
        {
            Field<GradType> neighbourSecondfGrad =
                refCast<const coupledFvPatchField<GradType> >
                (
                    secondfGrad.boundaryField()[patchI]
                ).patchNeighbourField();

            forAll(neighbourSecondfGrad, patchFaceI)
            {
                fGrad[FaceCells[patchFaceI]] +=
                    0.5*patchOwnerLs[patchFaceI]
                    *(
                        patchD[patchFaceI]
                      & (
                            neighbourSecondfGrad[patchFaceI]
                          - secondfGrad[FaceCells[patchFaceI]]
                        )
                    );
            }
        }
        else
        {
            const fvPatchField<GradType>& patchSecondfGrad =
                secondfGrad.boundaryField()[patchI];

            forAll(patchSecondfGrad, patchFaceI)
            {
                fGrad[FaceCells[patchFaceI]] +=
                    0.5*patchOwnerLs[patchFaceI]
                    *(
                        patchD[patchFaceI]
                      & (
                            patchSecondfGrad[patchFaceI]
                          - secondfGrad[FaceCells[patchFaceI]]
                        )
                    );
            }
        }
    }

    fGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, fGrad);

    return tfGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
