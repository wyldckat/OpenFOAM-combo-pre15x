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
    STARCoordinateRotation

Description
    The 3 rotations are defined in the STAR convention
    (around Z, around X' and around Y').


\*---------------------------------------------------------------------------*/

#ifndef STARCoordinateRotation_H
#define STARCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class STARCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class STARCoordinateRotation
:
    public coordinateRotation
{
    // Private data

        //- Local-to-global transformation tensor
        tensor R_;


public:

    // Constructors

        //- Construct from components
        STARCoordinateRotation
        (
            const scalar rotZ,
            const scalar rotY,
            const scalar rotX
        )
        :
            R_
            (
                cos(rotZ)*cos(rotY) + sin(rotZ)*sin(rotX)*sin(rotY),
                -sin(rotZ)*cos(rotX),
                -cos(rotZ)*sin(rotY) + sin(rotZ)*sin(rotX)*cos(rotY),

                sin(rotZ)*cos(rotY) - cos(rotZ)*sin(rotX)*sin(rotY),
                cos(rotZ)*cos(rotX),
                -sin(rotZ)*sin(rotY) - cos(rotZ)*sin(rotX)*cos(rotY),

                cos(rotX)*sin(rotY),
                sin(rotX),
                cos(rotX)*cos(rotY)
            )
        {}


    // Destructor

        virtual ~STARCoordinateRotation()
        {}


    // Return local-to-global transformation tensor
    virtual const tensor& R() const
    {
        return R_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
