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
    MUSCL

Description
    Class with weight function which returns the weighting factors for the
    van Leer's MUSCL differencing scheme.

    Used in conjunction with the template class NVDscheme.

SourceFiles
    MUSCLMake.C

\*---------------------------------------------------------------------------*/

#ifndef MUSCL_H
#define MUSCL_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MUSCLWeight Declaration
\*---------------------------------------------------------------------------*/

class MUSCLWeight
{

public:

    MUSCLWeight(Istream&)
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
        gradf = stabilise(gradf, SMALL);

        scalar r = 2*gradcf/gradf - 1;
        scalar limiter = max(min(min(2*r, 0.5*r + 0.5), 2), 0);

        return limiter*cdWeight + (1 - limiter)*udWeight;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
