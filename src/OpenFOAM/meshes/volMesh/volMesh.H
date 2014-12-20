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
    volMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

\*---------------------------------------------------------------------------*/

#ifndef volMesh_H
#define volMesh_H

#include "fvMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class volMesh Declaration
\*---------------------------------------------------------------------------*/

class volMesh
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

public:

    // Public typedefs

        typedef fvMesh Mesh;
        typedef fvBoundaryMesh BoundaryMesh;

    // Constructors

        //- Construct from fvMesh
        explicit volMesh(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Member Functions

        //- Return size
        label size() const
        {
            return size(mesh_);
        }

        //- Return size
        static label size(const Mesh& mesh)
        {
            return mesh.nCells();
        }

        //- Return cell centres
        const volVectorField& C()
        {
            return mesh_.C();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
