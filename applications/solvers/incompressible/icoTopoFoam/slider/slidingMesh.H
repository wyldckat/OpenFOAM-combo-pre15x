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
    slidingMesh

Description

SourceFiles
    slidingMeshI.H
    slidingMesh.C
    slidingMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef slidingMesh_H
#define slidingMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class slidingMesh Declaration
\*---------------------------------------------------------------------------*/

class slidingMesh
:
    public fvMesh
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        slidingMesh(const slidingMesh&);

        //- Disallow default bitwise assignment
        void operator=(const slidingMesh&);


public:

    // Constructors

        //- Construct from IOobject
        explicit slidingMesh(const IOobject& io)
        :
            fvMesh(io)
        {}


    // Destructor

        ~slidingMesh()
        {}


    // Member Functions

        //- Slide
        void slide()
        {
            vector motionVec(0.3, 0, 0);
//             vector motionVec(-0.3, 0, 0);

            movePoints
            (
                allPoints()
              + pos(allPoints().component(vector::Y) - 1.001)*motionVec
            );

            updateTopology();
            movePoints
            (
                oldAllPoints()
              + pos(oldAllPoints().component(vector::Y) - 1.001)*motionVec
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
