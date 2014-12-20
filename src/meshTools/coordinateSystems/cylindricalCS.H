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
    cylindricalCS

Description
    Cylindrical coordinate system: axis is the local z-axis, direction the
    local x-axis

SourceFiles
    cylindricalCS.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalCS_H
#define cylindricalCS_H

#include "cartesianCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cylindricalCS Declaration
\*---------------------------------------------------------------------------*/

class cylindricalCS
:
    public cartesianCS
{

public:

    //- Runtime type information
    TypeName("cylindrical");


    // Constructors

        //- Construct from components
        cylindricalCS
        (
            const word& name,
            const vector& origin,
            const vector& axis,
            const vector& direction
        );

        //- Construct from origin and rotation angles
        cylindricalCS
        (
            const word& name,
            const vector& origin,
            const coordinateRotation& cr
        );

        //- Construct from dictionary
        cylindricalCS(const word& name, const dictionary& dict);


    // Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        virtual vector toGlobal(const vector& localV) const;
        virtual tmp<vectorField> toGlobal(const vectorField& localV) const;

        //- Convert from global Cartesian coordinate system to the local system
        virtual vector toLocal(const vector& globalV) const;
        virtual tmp<vectorField> toLocal(const vectorField& globalV) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
