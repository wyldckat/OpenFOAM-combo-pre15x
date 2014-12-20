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
    gammaCompressionLimiter

Description
    Class with limiter function which returns the limiter for the
    TVD limited linear differencing scheme based on r obtained from the
    LimiterFunc class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    gammaCompression.C

\*---------------------------------------------------------------------------*/

#ifndef gammaCompression_H
#define gammaCompression_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gammaCompressionWeight Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class gammaCompressionLimiter
:
    public LimiterFunc
{

public:

    gammaCompressionLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        //return min(max(4*min(phiP*(1 - phiP), phiN*(1 - phiN)), 0), 1);
        return min(max(1 - max(sqr(1 - 4*phiP*(1 - phiP)), sqr(1 - 4*phiN*(1 - phiN))), 0), 1);
        //return min(max(1 - max(sqr(sqr(1 - 4*phiP*(1 - phiP))), sqr(sqr(1 - 4*phiN*(1 - phiN)))), 0), 1);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
