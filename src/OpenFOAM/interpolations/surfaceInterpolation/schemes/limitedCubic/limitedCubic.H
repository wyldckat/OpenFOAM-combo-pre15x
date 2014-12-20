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
    limitedCubic

Description
    Class with weight function which returns the weighting factors for the
    TVD limited centred-cubic differencing scheme.

    Used in conjunction with the template class NVDscheme.

SourceFiles
    limitedCubicMake.C

\*---------------------------------------------------------------------------*/

#ifndef limitedCubic_H
#define limitedCubic_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class limitedCubicWeight Declaration
\*---------------------------------------------------------------------------*/

class limitedCubicWeight
{
    scalar kl_, ku_;

public:

    limitedCubicWeight(Istream&)
    {}

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
        scalar phiU;

        if (faceFlux > 0)
        {
            gradcf = dHat & gradcP;
            udWeight = 1;
            phiU = phiP;
        }
        else
        {
            gradcf = dHat & gradcN;
            udWeight = 0;
            phiU = phiN;
        }

        // Stabilise for division
        gradf = stabilise(gradf, SMALL);

        scalar r = 2*gradcf/gradf - 1;

        // Calculate the face value using cubic interpolation
        scalar phif = 
            cdWeight*(phiP - 0.25*(d & gradcN))
          + (1 - cdWeight)*(phiN + 0.25*(d & gradcP));

        scalar phiCD = cdWeight*phiP + (1 - cdWeight)*phiN;

        // Calculate the effective limiter for the cubic interpolation
        scalar cubicLimiter = (phif - phiU)/stabilise(phiCD - phiU, SMALL);

        // ad-hoc limit on the TVD parameter to only limit the weighting
        // factors when absolutly necessary
        r = max(r, 0);
        if (r > 0.05) r = 1000;

        // Limit the limiter to obey the TVD constraint
        scalar limiter = max(min(min(2*r, cubicLimiter), 2), 0);

        //scalar limiter =
        //    0.5*max(min(min(2*r, cubicLimiter), 2), 0)
        //  + 0.5*max(min(cubicLimiter, 2);

        return limiter*cdWeight + (1 - limiter)*udWeight;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
