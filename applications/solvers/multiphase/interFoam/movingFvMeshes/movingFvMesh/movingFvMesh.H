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
    movingFvMesh

Description

SourceFiles
    movingFvMesh.C
    newMovingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef movingFvMesh_H
#define movingFvMesh_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class movingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class movingFvMesh
:
    public fvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingFvMesh(const movingFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const movingFvMesh&);


public:

    //- Runtime type information
    TypeName("movingFvMesh");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            movingFvMesh,
            IOobject,
            (const IOobject& io),
            (io)
        );


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit movingFvMesh(const IOobject& io);


    // Selectors

        //- Select null constructed
        static autoPtr<movingFvMesh> New(const IOobject& io);


    // Destructor

        virtual ~movingFvMesh();


    // Member Functions

        // Edit

            virtual void move() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
