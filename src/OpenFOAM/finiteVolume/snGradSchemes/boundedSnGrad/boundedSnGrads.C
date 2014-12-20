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
    Central-difference snGrad scheme with bounded non-orthogonal correction.
    Boundedness is handled by converting the non-orthogonal part into a
    convection term and then applying upwind-differencing to it

\*---------------------------------------------------------------------------*/

#include "boundedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "upwind.H"
#include "correctedSnGrad.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeSnGradTypeScheme(boundedSnGrad, scalar)


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundedSnGrad<scalar>::~boundedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<GeometricField<scalar, fvPatchField, surfaceMesh> >
boundedSnGrad<scalar>::correction
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    surfaceScalarField upwindValue = linearInterpolate(fvc::average(vf));

    // stabilise the denominator
    upwindValue +=
        sign(upwindValue)
       *dimensionedScalar("small", upwindValue.dimensions(), SMALL);

    // calculate cross-diffusion flux
    surfaceScalarField crossFlux = 
        correctedSnGrad<scalar>(mesh()).correction(vf)/upwindValue;

    // Upwind vf according to the cross-diffusion flux
    // note: -crossFlux is used for upwinding because the
    //   laplacian matrix is negated when it is added to the convection
    //   term.  If crossFlux were used the term would be down-winded
    return crossFlux*upwind<scalar>(vf.mesh(), -crossFlux).interpolate(vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
