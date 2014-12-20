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
    Phi

Description
    Class with weight function which returns the weighting factors for the
    standard Phi differencing scheme.

    Used in conjunction with the template class PhiScheme.

SourceFiles
    Phi.C

\*---------------------------------------------------------------------------*/

#ifndef Phi_H
#define Phi_H

#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PhiWeight Declaration
\*---------------------------------------------------------------------------*/

class PhiWeight
{
    scalar k_;

public:

    PhiWeight(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("PhiWeight(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }
    }


    scalar weight
    (
        scalar cdWeight,
        scalar faceFlux,
        const vector& phiP,
        const vector& phiN,
        const vector& Sf,
        const scalar&
    ) const
    {
        scalar weight = 
            k_*((faceFlux - (Sf&phiN))/stabilise(Sf&(phiP - phiN), SMALL))
          + (1.0 - k_)*cdWeight;

        if (faceFlux > 0)
        {
            weight = min(max(weight, cdWeight), 1);
        }
        else
        {
            weight = min(max(weight, 0), cdWeight);
        }

        return weight;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
