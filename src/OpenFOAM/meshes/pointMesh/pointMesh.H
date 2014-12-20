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
    pointMesh

Description

\*---------------------------------------------------------------------------*/

#ifndef pointMesh_H
#define pointMesh_H

#include "polyMesh.H"
#include "pointBoundaryMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pointMesh Declaration
\*---------------------------------------------------------------------------*/

class pointMesh
{
    // Permanent data

        //- Reference to underlying polyMesh
        const polyMesh& mesh_;

        //- Boundary mesh
        pointBoundaryMesh boundary_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointMesh(const pointMesh&);

        //- Disallow default bitwise assignment
        void operator=(const pointMesh&);

        //- Check for and create a parallel point patch
        void addParallelPointPatch();


public:

    typedef pointMesh Mesh;
    typedef pointBoundaryMesh BoundaryMesh;

    // Constructors

        //- Construct from polyMesh
        explicit pointMesh(const polyMesh& pMesh);


    // Member Functions

        //- Return reference to polyMesh
        const polyMesh& operator()() const
        {
            return mesh_;
        }

        //- Return number of points
        label size() const
        {
            return size(*this);
        }

        //- Return number of points
        static label size(const Mesh& mesh)
        {
            return mesh.mesh_.nPoints();
        }

        //- Return number of points (checking of patch size.  Reconsider)
        label nPoints() const
        {
            return mesh_.nPoints();
        }

        //- Return number of cells (checking of patch size.  Reconsider)
        label nCells() const
        {
            return mesh_.nCells();
        }

        //- Return reference to boundary mesh
        const pointBoundaryMesh& boundary() const
        {
            return boundary_;
        }

        //- Correct stuff for moving mesh.
        bool movePoints()
        {
            return true;
        }

        //- Return parallel info
        const parallelInfo& parallelData() const
        {
            return mesh_.parallelData();
        }


    // Member Operators

        bool operator!=(const pointMesh& pm) const
        {
            return &pm != this;
        }

        bool operator==(const pointMesh& pm) const
        {
            return &pm == this;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
