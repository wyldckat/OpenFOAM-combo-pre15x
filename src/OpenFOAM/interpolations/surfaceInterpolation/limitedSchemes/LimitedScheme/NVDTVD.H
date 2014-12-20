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

Class
    NVDTVD

Description

SourceFiles
    NVDTVD.C

\*---------------------------------------------------------------------------*/

#ifndef NVDTVD_H
#define NVDTVD_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LimitedScheme Declaration
\*---------------------------------------------------------------------------*/

class NVDTVD
{

public:

    typedef scalar phiType;
    typedef vector gradPhiType;


    // Null Constructor

        NVDTVD()
        {}


    // Member Functions

        scalar phict
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const scalar phiP,
            const scalar phiN,
            const vector& gradcP,
            const vector& gradcN,
            const vector& d
        ) const
        {
            scalar magd = mag(d);
            vector dHat = d/magd;

            scalar gradf = (phiN - phiP)/magd;

            scalar gradcf;

            if (faceFlux > 0)
            {
                gradcf = dHat & gradcP;
            }
            else
            {
                gradcf = dHat & gradcN;
            }

            // Stabilise for division
            gradcf = stabilise(gradcf, SMALL);

            return 1 - 0.5*gradf/gradcf;
        }


        scalar r
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const scalar phiP,
            const scalar phiN,
            const vector& gradcP,
            const vector& gradcN,
            const vector& d
        ) const
        {
            scalar magd = mag(d);
            vector dHat = d/magd;
        
            scalar gradf = (phiN - phiP)/magd;

            scalar gradcf;

            if (faceFlux > 0)
            {
                gradcf = dHat & gradcP;
            }
            else
            {
                gradcf = dHat & gradcN;
            }

            // Stabilise for division
            gradf = stabilise(gradf, SMALL);

            return 2*gradcf/gradf - 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
