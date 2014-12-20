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
    cubic

Description
    Class with operator() which returns the weighting factors for the
    cubic-centred differencing scheme.  Note that the weighting factors are
    not bounded between upwind and central-differencing, some downwind
    contribution is possible although the interpolate is limited to be between
    the upwind and downwind cell values.

    Used in conjunction with the template class NVDscheme although this scheme
    is not NVD.

SourceFiles
    cubicMake.C

\*---------------------------------------------------------------------------*/

#ifndef cubic_H
#define cubic_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cubicWeight Declaration
\*---------------------------------------------------------------------------*/

class cubicWeight
{

public:

    cubicWeight(Istream&)
    {}

    scalar weight
    (
        scalar cdWeight,
        scalar,
        scalar phiP,
        scalar phiN,
        const vector& gradcP,
        const vector& gradcN,
        const vector& d
    ) const
    {
        scalar phif = 
            cdWeight*(phiP - 0.25*(d & gradcN))
          + (1 - cdWeight)*(phiN + 0.25*(d & gradcP));

        // Limit the estimated face value between the upwind and downwind cell
        // values
        phif = min(phif, max(phiN, phiP));
        phif = max(phif, min(phiN, phiP));

        return (phif - phiN)/stabilise(phiP - phiN, SMALL);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
