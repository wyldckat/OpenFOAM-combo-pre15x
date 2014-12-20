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

Description
    Function to calculate the laminar flame speed of iso-octane,
    propane or methanol as a function of pressure, temperature and
    equivalence ratio.

    The correlation functions used are those proposed by Metgalchi
    and Keck in combustion and flame vol. 48, pp 191-210.

    The correlations are valid for
        pressure p = 0.4 - 50 atm
        unburnt gas temperature tu = 300 - 700 k
        equivalence ratio eqfur = 0.8 - 1.5

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "Sl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar Sl(scalar p, scalar t, scalar eqfur, fuelType iftype)
{
    scalar eqfurm = 1.0, bm = 0.0, b2 = 0.0;

    // convert total pressure from pa into atm
    scalar ptot = p/1.013e5;

    switch(iftype)
    {
        case ISOOCTANE:
            eqfurm = 1.13;
            bm = 0.2632;
            b2 = -0.8472;
        break;

        case PROPANE:
            eqfurm = 1.08;
            bm = 0.3422;
            b2 = -1.3865;
        break;

        case METHANOL:
            eqfurm = 1.11;
            bm = 0.3692;
            b2 = -1.4051;
        break;

        default:

            FatalErrorIn("Sl(scalar p, scalar t, eqfur, fuelType iftype)")
                << "unknown fuel type"
                << abort(FatalError);
    }


    scalar lflsp0 = bm + b2*(eqfur - eqfurm)*(eqfur - eqfurm);
    scalar var1 = (eqfur - 1.0);
    scalar alfa = 2.18 - 0.8*var1;
    scalar betaSl = -0.16 + 0.22*var1;

    return lflsp0 * pow(t/298.0, alfa)*pow(ptot, betaSl);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
