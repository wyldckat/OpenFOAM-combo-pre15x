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
    morphMesh

Description
    Wrapper around polyMesh to make public all the nice protected morphing
    functionality.

SourceFiles
    morphMesh.C

\*---------------------------------------------------------------------------*/

#ifndef morphMesh_H
#define morphMesh_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class Time;
class polyTopoChange;

/*---------------------------------------------------------------------------*\
                           Class morphMesh Declaration
\*---------------------------------------------------------------------------*/

class morphMesh
:
    public polyMesh
{

public:

    // Constructors

        //- Construct from IOobject
        explicit morphMesh(const IOobject& io);


        //- Construct from components without boundary.
        // Boundary is added using addPatches() member function
        morphMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );


    // Destructor

        ~morphMesh();


    // Member Functions

        //- Set the instance for mesh files
        void setInstance(const fileName& inst)
        {
            polyMesh::setInstance(inst);
        }

        //- Remove existing patches
        void removeBoundary()
        {
            polyMesh::removeBoundary();
        }

        //- Add boundary patches. Constructor helper
        void addPatches(const List<polyPatch*>& patches)
        {
            polyMesh::addPatches(patches);
        }

        //- Return raw points
        const pointField& allPoints() const
        {
            return polyMesh::allPoints();
        }

        //- Return raw faces
        const faceList& allFaces() const
        {
            return polyMesh::allFaces();
        }

        //- Return raw cells
        const cellList& allCells() const
        {
            return polyMesh::allCells();
        }

        //- Return old mesh motion points
        const pointField& oldAllPoints() const
        {
            return polyMesh::oldAllPoints();
        }

        //- Return face owner
        const labelList& allOwner() const
        {
            return polyMesh::allOwner();
        }

        //- Return face neighbour
        const labelList& allNeighbour() const
        {
            return polyMesh::allNeighbour();
        }

        //- Change topology of the mesh given a topological change
        virtual void updateTopology(const polyTopoChange& changer)
        {
            polyMesh::updateTopology(changer);
        }

        virtual bool setMorphTimeIndex(const label newTimeIndex) const
        {
            return polyMesh::setMorphTimeIndex(newTimeIndex);
        }

        //- Inplace modify mesh primitives. Use addPatches to add boundary
        //  later.
        void reset
        (
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        )
        {
            polyMesh::reset(points, faces, cells);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
