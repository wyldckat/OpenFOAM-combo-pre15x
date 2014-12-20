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
    
\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "limitedGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvGradScheme(limitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Limited scalar gradient
template<>
tmp<volVectorField> limitedGrad<scalar>::grad
(
    const volScalarField& vsf
)
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volVectorField> tGrad = basicGradScheme_().grad(vsf);

    // Limit the gradient

    // Owner/neighbour addressing
    volVectorField& g = tGrad();

    const unallocLabelList& Owner = mesh.owner();
    const unallocLabelList& Neighbour = mesh.neighbour();

    const vectorField& cellCentres = vsf.mesh().C();

    const vectorField& faceCentres = vsf.mesh().Cf().internalField();

    // create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    forAll(Owner, faceI)
    {
        label own = Owner[faceI];

        label nei = Neighbour[faceI];

        scalar vsfOwn = vsf[own];

        scalar vsfNei = vsf[nei];

        scalar maxFace = max(vsfOwn, vsfNei);

        scalar minFace = min(vsfOwn, vsfNei);

        // Owner side
        scalar extrapolatedOwner =
            (faceCentres[faceI] - cellCentres[own]) & g[own];

        if (vsfOwn + extrapolatedOwner > maxFace + SMALL)
        {
            limiter[own] =
                min
                (
                    limiter[own],
                    (maxFace - vsfOwn)/extrapolatedOwner
                );
        }
        else if (vsfOwn + extrapolatedOwner < minFace - SMALL)
        {
            limiter[own] =
                min
                (
                    limiter[own],
                    (minFace - vsfOwn)/extrapolatedOwner
                );
        }

        // neighbour side
        scalar extrapolatedNeighbour =
            (faceCentres[faceI] - cellCentres[nei]) & g[nei];

        if (vsfNei + extrapolatedNeighbour > maxFace + SMALL)
        {
            limiter[nei] =
                min
                (
                    limiter[nei],
                     (maxFace - vsfNei)/extrapolatedNeighbour
                );
        }
        else if (vsfNei + extrapolatedNeighbour < minFace - SMALL)
        {
            limiter[nei] =
                min
                (
                    limiter[nei],
                    (minFace - vsfNei)/extrapolatedNeighbour
                );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << max(limiter)
            << " min = " << min(limiter)
            << " average: " << Foam::average(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
tmp<volTensorField> limitedGrad<vector>::grad
(
    const volVectorField& vvf
)
{
    const fvMesh& mesh = vvf.mesh();

    tmp<volTensorField> tGrad = basicGradScheme_().grad(vvf);

    // Limit the gradient

    // Owner/neighbour addressing
    volTensorField& g = tGrad();

    const unallocLabelList& Owner = mesh.owner();
    const unallocLabelList& Neighbour = mesh.neighbour();

    const vectorField& cellCentres = vvf.mesh().C().internalField();

    const vectorField& faceCentres = vvf.mesh().Cf().internalField();

    // create limiter
    scalarField limiter(vvf.internalField().size(), 1.0);

    forAll(Owner, faceI)
    {
        label own = Owner[faceI];
        label nei = Neighbour[faceI];

        vector vvfOwn = vvf[own];
        vector vvfNei = vvf[nei];

        vector gradfHat = (vvfNei - vvfOwn);
        gradfHat /= mag(gradfHat) + SMALL;

        scalar vsfOwn = gradfHat & vvfOwn;
        scalar vsfNei = gradfHat & vvfNei;

        scalar maxFace = max(vsfOwn, vsfNei);
        scalar minFace = min(vsfOwn, vsfNei);

        // Owner side
        scalar extrapolatedOwner =
            gradfHat & ((faceCentres[faceI] - cellCentres[own]) & g[own]);

        if (vsfOwn + extrapolatedOwner > maxFace + SMALL)
        {
            limiter[own] =
                min
                (
                    limiter[own],
                    (maxFace - vsfOwn)/extrapolatedOwner
                );
        }
        else if (vsfOwn + extrapolatedOwner < minFace - SMALL)
        {
            limiter[own] =
                min
                (
                    limiter[own],
                    (minFace - vsfOwn)/extrapolatedOwner
                );
        }


        // neighbour side
        scalar extrapolatedNeighbour =
            gradfHat & ((faceCentres[faceI] - cellCentres[nei]) & g[nei]);

        if (vsfNei + extrapolatedNeighbour > maxFace + SMALL)
        {
            limiter[nei] =
                min
                (
                    limiter[nei],
                    (maxFace - vsfNei)/extrapolatedNeighbour
                );
        }
        else if (vsfNei + extrapolatedNeighbour < minFace - SMALL)
        {
            limiter[nei] =
                min
                (
                    limiter[nei],
                    (minFace - vsfNei)/extrapolatedNeighbour
                );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vvf.name()
            << " max = " << max(limiter)
            << " min = " << min(limiter)
            << " average: " << Foam::average(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vvf, g);

    return tGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
