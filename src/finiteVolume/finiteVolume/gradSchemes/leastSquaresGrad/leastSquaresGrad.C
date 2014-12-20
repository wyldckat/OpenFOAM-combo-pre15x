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

#include "leastSquaresGrad.H"
#include "leastSquaresVectors.H"
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
) const
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
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);
    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    // Owner/neighbour addressing
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    forAll(own, faceI)
    {
        register label ownFaceI = own[faceI];
        register label neiFaceI = nei[faceI];

        Type deltaVsf = vsf[neiFaceI] - vsf[ownFaceI];

        lsGrad[ownFaceI] += ownLs[faceI]*deltaVsf;
        lsGrad[neiFaceI] -= neiLs[faceI]*deltaVsf;
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchI)
    {
        const fvPatchVectorField& patchOwnLs =
            ownLs.boundaryField()[patchI];

        const labelList::subList FaceCells =
            lsGrad.boundaryField()[patchI].patch().faceCells();

        if (vsf.boundaryField()[patchI].coupled())
        {
            Field<Type> neiVsf =
                refCast<const coupledFvPatchField<Type> >
                (
                    vsf.boundaryField()[patchI]
                ).patchNeighbourField();

            forAll(neiVsf, patchFaceI)
            {
                lsGrad[FaceCells[patchFaceI]] +=
                    patchOwnLs[patchFaceI]
                   *(neiVsf[patchFaceI] - vsf[FaceCells[patchFaceI]]);
            }
        }
        else
        {
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchI];

            forAll(patchVsf, patchFaceI)
            {
                lsGrad[FaceCells[patchFaceI]] +=
                     patchOwnLs[patchFaceI]
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
