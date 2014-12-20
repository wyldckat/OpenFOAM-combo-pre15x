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
    Central-difference snGrad scheme with limited non-orthogonal correction.
    The limiter is controlled by a coefficient with a value between 0 and 1
    which when 0 switches the correction off and the scheme behaves as
    uncorrectedSnGrad, when set to 1 the full correction is applied and the
    scheme behaves as correctedSnGrad and when set to 0.5 the limiter is
    calculated such that the non-orthogonal contribution does not exceed the
    orthogonal part.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "limitedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "correctedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
limitedSnGrad<Type>::~limitedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
limitedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // The limiter on the non-orthogonal correction is such that the
    // non-orthogonal correction at most equals the orthogonal term
    // (its direction cannot be checked in for anything else than a scalar)

    GeometricField<Type, fvPatchField, surfaceMesh> corr = 
        correctedSnGrad<Type>(this->mesh()).correction(vf);

    surfaceScalarField limiter
    (
        min
        (
            (
                limitCoeff_*
                mag
                (
                    snGradScheme<Type>::snGrad(vf, deltaCoeffs(vf))
                )
              + dimensionedScalar("small", corr.dimensions(), SMALL)
            )
           /(
                (1 - limitCoeff_)*mag(corr)
              + dimensionedScalar("small", corr.dimensions(), SMALL)
            ),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // for faces where the delta coefficient is negative
    // (non-orthoginality over 90 deg), set the limiter to zero
    limiter.internalField() *= pos(deltaCoeffs(vf)().internalField());

    if (fv::debug)
    {
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
    }

    return limiter*corr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
