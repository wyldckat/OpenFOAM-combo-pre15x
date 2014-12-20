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
    subCycle

Description

\*---------------------------------------------------------------------------*/

#ifndef subCycle_H
#define subCycle_H

#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class subCycle Declaration
\*---------------------------------------------------------------------------*/

template<class GeometricField>
class subCycle
:
    public subCycleTime
{
    // Private data

        GeometricField& gf_;
        //GeometricField  gf0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        subCycle(const subCycle<GeometricField>&);

        //- Disallow default bitwise assignment
        void operator=(const subCycle<GeometricField>&);


public:

    // Constructors

        //- Construct field and number of sub-cycles
        subCycle(GeometricField& gf, const label nSubCycles)
        :
            subCycleTime(const_cast<Time&>(gf.time()), nSubCycles),
            gf_(gf)
            //, gf0_(gf.oldTime())
        {
            gf.oldTime();
        }


    // Destructor

        ~subCycle()
        {
            //gf_.oldTime() = gf0_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
