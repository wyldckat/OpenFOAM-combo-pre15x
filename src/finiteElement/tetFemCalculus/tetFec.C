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
    Class of static functions to calculate explicit finite element derivatives.

\*---------------------------------------------------------------------------*/

#include "tetCell.H"
#include "tetCellList.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        tetPolyPatchField,
        tetPointMesh
    >
> 
tetFec::grad
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& psi
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const tetPolyMesh& tetMesh = psi.mesh();

    tmp<GeometricField<GradType, tetPolyPatchField, tetPointMesh> > tFemGrad
    (
        new GeometricField<GradType, tetPolyPatchField, tetPointMesh>
        (
            IOobject
            (
                "grad("+psi.name()+')',
                psi.instance(),
                psi.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh,
            dimensioned<GradType>
            (
                "zero",
                psi.dimensions()/dimLength,
                pTraits<GradType>::zero
            )
        )
    );

    GeometricField<GradType, tetPolyPatchField, tetPointMesh>& femGrad = 
        tFemGrad();

    pointField points = tetMesh.points();
    cellShapeList tetCellShapes = tetMesh.tetCells();

    scalarField weights (points.size(), 0.0);

    forAll (tetCellShapes, tetI)
    {
        cellShape& curShape = tetCellShapes[tetI];

        tetPointRef curTetrahedron
        (
            points[curShape[0]],
            points[curShape[1]],
            points[curShape[2]],
            points[curShape[3]]
        );

        GradType tetGrad = 
          - (1.0/3.0)*
            (
                curTetrahedron.Sa()*psi.internalField()[curShape[0]]
              + curTetrahedron.Sb()*psi.internalField()[curShape[1]]
              + curTetrahedron.Sc()*psi.internalField()[curShape[2]]
              + curTetrahedron.Sd()*psi.internalField()[curShape[3]]
            )/curTetrahedron.mag();

        forAll (curShape, pointI)
        {
            scalar weight = 
                curTetrahedron.mag()/
                mag
                (
                    points[curShape[pointI]]
                  - curShape.centre(points)
                );

            femGrad.internalField()[curShape[pointI]] += weight*tetGrad;

            weights[curShape[pointI]] += weight;
        }
    }

    femGrad.internalField() /= weights;

    return tFemGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
