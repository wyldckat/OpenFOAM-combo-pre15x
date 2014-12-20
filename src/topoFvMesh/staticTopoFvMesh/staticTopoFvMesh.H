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
    staticTopoFvMesh

Description

SourceFiles
    staticTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef staticTopoFvMesh_H
#define staticTopoFvMesh_H

#include "topoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class staticTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class staticTopoFvMesh
:
    public topoFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        staticTopoFvMesh(const staticTopoFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const staticTopoFvMesh&);


public:

    //- Runtime type information
    TypeName("staticTopoFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit staticTopoFvMesh(const IOobject& io);


    // Destructor

        ~staticTopoFvMesh();


    // Member Functions

        //- Move and morph
        virtual void moveAndMorph();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
