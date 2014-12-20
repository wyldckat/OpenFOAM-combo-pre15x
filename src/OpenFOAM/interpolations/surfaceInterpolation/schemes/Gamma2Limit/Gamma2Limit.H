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

Class
    Gamma2Limit

Description
    Class with operator() which returns the weighting factors for the
    Gamma2Limit differencing scheme.  Used in conjunction with the template
    class NVDscheme.  The scheme is a Gamma2 scheme with a defined lower and
    upper bound for the variable that is being solved.  If any of the
    bounds are violated, the scheme recovers to upwind in order to
    recover boundedness.

SourceFiles
    Gamma2Limit.C

\*---------------------------------------------------------------------------*/

#ifndef Gamma2Limit_H
#define Gamma2Limit_H

#include "scalar.H"
#include "vector.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Gamma2LimitWeight Declaration
\*---------------------------------------------------------------------------*/

class Gamma2LimitWeight
{
private:

    //- Lower and upper bound of the variable
    scalar lowerBound_, upperBound_;

    //- Control factor of the Gamma2 scheme
    scalar k_;

    void checkParameters(Istream& is)
    {
        if (lowerBound_ > upperBound_)
        {
            FatalIOErrorIn("checkParameters()", is)
                << "Invalid bounds.  Lower = " << lowerBound_
                << "  Upper = " << upperBound_
                << ".  Lower bound is higher than the upper bound."
                << exit(FatalIOError);
        }

        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("checkParameters()", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        // Rescale k_ to be >= 0 and <= 0.5 (TVD conformant)
        // and avoid the /0 when k_ = 0
        k_ = max(k_/2.0, SMALL);
    }


public:

    Gamma2LimitWeight
    (
        const scalar lowerBound,
        const scalar upperBound,
        Istream& is
    )
    :
        lowerBound_(lowerBound),
        upperBound_(upperBound),
        k_(readScalar(is))
    {
        checkParameters(is);
    }

    Gamma2LimitWeight(Istream& is)
    :
        lowerBound_(readScalar(is)),
        upperBound_(readScalar(is)),
        k_(readScalar(is))
    {
        checkParameters(is);
    }


    scalar weight
    (
        scalar cdWeight,
        scalar faceFlux,
        scalar phiP,
        scalar phiN,
        const vector& gradcP,
        const vector& gradcN,
        const vector& d
    ) const
    {
        scalar magd = mag(d);
        vector dHat = d/mag(d);

        scalar gradf = (phiN - phiP)/magd;

        scalar gradcf;
        scalar udWeight;

        if (faceFlux > 0)
        {
            gradcf = dHat & gradcP;
            udWeight = 1;
        }
        else
        {
            gradcf = dHat & gradcN;
            udWeight = 0;
        }

        // Stabilise for division
        gradcf = stabilise(gradcf, SMALL);

        scalar phict = 1 - 0.5*gradf/gradcf;
        scalar limiter = 1 - sqr(1 - min(max(phict/k_, 0), 1));

        // If the field values are not between 0 and 1 use upwind
        if
        (
            phiP < lowerBound_ || phiP > upperBound_
         || phiN < lowerBound_ || phiN > upperBound_
        )
        {
            limiter = 0;
        }

        return limiter*cdWeight + (1 - limiter)*udWeight;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
