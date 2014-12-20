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

#include "leastSquaresGrad.H"
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
leastSquaresGrad<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
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
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad();

    // Get reference to least square vectors
    const surfaceVectorField& ownerLs =
        mesh.surfaceInterpolation::leastSquarePvectors();

    const surfaceVectorField& neighbourLs =
        mesh.surfaceInterpolation::leastSquareNvectors();

    // Owner/neighbour addressing
    const unallocLabelList& Owner = mesh.owner();
    const unallocLabelList& Neighbour = mesh.neighbour();

    forAll(Owner, faceI)
    {
        lsGrad[Owner[faceI]] +=
            ownerLs[faceI]*(vsf[Neighbour[faceI]] - vsf[Owner[faceI]]);

        lsGrad[Neighbour[faceI]] -=
            neighbourLs[faceI]*(vsf[Neighbour[faceI]] - vsf[Owner[faceI]]);
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchI)
    {
        const fvPatchVectorField& patchOwnerLs =
            ownerLs.boundaryField()[patchI];

        const labelList::subList FaceCells =
            lsGrad.boundaryField()[patchI].patchMesh().faceCells();

        if (vsf.boundaryField()[patchI].coupled())
        {
            Field<Type> neighbourVsf =
                refCast<const coupledFvPatchField<Type> >
                (
                    vsf.boundaryField()[patchI]
                ).patchNeighbourField();

            forAll(neighbourVsf, patchFaceI)
            {
                lsGrad[FaceCells[patchFaceI]] +=
                    patchOwnerLs[patchFaceI]
                    *(neighbourVsf[patchFaceI] - vsf[FaceCells[patchFaceI]]);
            }
        }
        else
        {
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchI];

            forAll(patchVsf, patchFaceI)
            {
                lsGrad[FaceCells[patchFaceI]] +=
                    patchOwnerLs[patchFaceI]
                    *(patchVsf[patchFaceI] - vsf[FaceCells[patchFaceI]]);
            }
        }
    }

    lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
