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
    movingInkJetFvMesh

Description

SourceFiles
    movingInkJetFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef movingInkJetFvMesh_H
#define movingInkJetFvMesh_H

#include "movingFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class movingInkJetFvMesh Declaration
\*---------------------------------------------------------------------------*/

class movingInkJetFvMesh
:
    public movingFvMesh
{
    // Private data

        dictionary movingMeshCoeffs_;

        scalar amplitude_;
        scalar frequency_;
        scalar refPlaneX_;

        pointIOField stationaryPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingInkJetFvMesh(const movingInkJetFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const movingInkJetFvMesh&);


public:

    //- Runtime type information
    TypeName("movingInkJetFvMesh");


    // Constructors

        //- Construct from IOobject
        movingInkJetFvMesh(const IOobject& io);


    // Destructor

        ~movingInkJetFvMesh();


    // Member Functions

        // Edit

            void move();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
