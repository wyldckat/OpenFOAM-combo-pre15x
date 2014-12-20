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

#include "multivariateSelectionScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void multivariateSelectionScheme<Type>::makeWeights()
{
    const fvMesh& mesh = faceFlux_.mesh();

    const surfaceScalarField& linearWeights = mesh.surfaceInterpolation::weights();

    typename multivariateSurfaceInterpolationScheme<Type>::
        fieldTable::const_iterator iter = this->fields().begin();

    weights_ = upwind<Type>(mesh, faceFlux_).weights(*iter());
    weights_.boundaryField() += SMALL;

    surfaceScalarField limiter =
    (
        weights_
      - surfaceInterpolationScheme<Type>::New
        (
            mesh,
            faceFlux_,
            schemes_.lookup(iter()->name())
        )().weights(*iter())
    )/(weights_ - linearWeights);

    for (++iter; iter != this->fields().end(); ++iter)
    {
        limiter = min
        (
            limiter,
            (
                weights_
              - surfaceInterpolationScheme<Type>::New
                (
                    mesh,
                    faceFlux_,
                    schemes_.lookup(iter()->name())
                )().weights(*iter())
            )/(weights_ - linearWeights)
        );
    }

    weights_ = limiter*linearWeights + (1.0 - limiter)*weights_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
