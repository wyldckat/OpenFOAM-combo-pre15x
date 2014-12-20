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
    Gamma2V

Description
    Class with weight function which returns the weighting factors for the
    Gamma2V differencing scheme.  Used in conjunction with the template class
    NVDscheme.

SourceFiles
    Gamma2V.C

\*---------------------------------------------------------------------------*/

#ifndef Gamma2V_H
#define Gamma2V_H

#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Gamma2VWeight Declaration
\*---------------------------------------------------------------------------*/

class Gamma2VWeight
{
    scalar k_;

public:

    Gamma2VWeight(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("Gamma2VWeight(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        // Rescale k_ to be >= 0 and <= 0.5 (TVD conformant)
        // and avoid the /0 when k_ = 0
        k_ = max(k_/2.0, SMALL);
    }


    scalar weight
    (
        scalar cdWeight,
        scalar faceFlux,
        const vector& phiP,
        const vector& phiN,
        const tensor& gradcP,
        const tensor& gradcN,
        const vector& d
    ) const
    {
        scalar magd = mag(d);
        vector dHat = d/mag(d);

        vector gradfV = (phiN - phiP)/magd;
        scalar gradf = mag(gradfV);

        vector gradfHat = gradfV/(gradf + SMALL);

        scalar gradcf;
        scalar udWeight;

        if (faceFlux > 0)
        {
            gradcf = gradfHat & (dHat & gradcP);
            udWeight = 1;
        }
        else
        {
            gradcf = gradfHat & (dHat & gradcN);
            udWeight = 0;
        }

        // Stabilise for division
        gradcf = stabilise(gradcf, SMALL);

        scalar phict = 1 - 0.5*gradf/gradcf;
        scalar limiter = 1 - sqr(1 - min(max(phict/k_, 0), 1));

        return limiter*cdWeight + (1 - limiter)*udWeight;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
