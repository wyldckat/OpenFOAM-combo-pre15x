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
    Sl

Description
    Function to calculate the laminar flame speed of iso-octane,
    propane or methanol as a function of pressure, temperature and
    equivalence ratio.  The correlations used are those of Keck.

    The correlation functions used are those proposed by Metgalchi
    and Keck in combustion and flame vol. 48, pp 191-210.

    The correlations are valid for
        pressure p = 0.4 - 50 atm
        unburnt gas temperature tu = 300 - 700 k
        equivalence ratio eqfur = 0.8 - 1.5

SourceFiles
    Sl.C
    SlIO.C


\*---------------------------------------------------------------------------*/

#ifndef Sl_H
#define Sl_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// enumeration definitions

enum fuelType
{
    ISOOCTANE,
    PROPANE,
    METHANOL
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//scalar Sl(scalar p, scalar t, scalar eqfur, fuelType iftype);

#include "Sl.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
