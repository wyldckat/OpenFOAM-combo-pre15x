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
    
\*---------------------------------------------------------------------------*/

#include "gaussGrad.H"
#include "fvcGrad.H"
#include "GeometricField.H"

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
gaussGrad<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        fvc::surfaceIntegrate
        (
            vsf.mesh().Sf()
           *tinterpScheme_().interpolate(vsf)
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    gGrad.rename("grad(" + vsf.name() + ')');
    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    // Henry's boundary gradient trick
    forAll(vsf.boundaryField(), patchI)
    {
        if (!vsf.boundaryField()[patchI].coupled())
        {
            vectorField n =
                vsf.mesh().Sf().boundaryField()[patchI]
               /vsf.mesh().magSf().boundaryField()[patchI];

            gGrad.boundaryField()[patchI] += n *
            (
                vsf.boundaryField()[patchI].snGrad()
              - (n & gGrad.boundaryField()[patchI])
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
