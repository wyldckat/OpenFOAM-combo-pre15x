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
    morphFvMesh

Description
    Makes topo change functionality of fvMesh public.

SourceFiles
    morphFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef morphFvMesh_H
#define morphFvMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyTopoChange;
class IOobjectList;

/*---------------------------------------------------------------------------*\
                           Class morphFvMesh Declaration
\*---------------------------------------------------------------------------*/

class morphFvMesh
:
    public fvMesh
{

    //- Helper function to read fields
    template<class GeoField>
    void readFields(const IOobjectList&, ptrList<GeoField>&) const;

public:

    // Constructors

        //- Construct from IOobject
        explicit morphFvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        morphFvMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );


    // Destructor

        ~morphFvMesh();


    // Member Functions

        //- Add boundary patches. Constructor helper
        void addFvPatches(const List<polyPatch*>& patches);

        //- Remove fvPatches
        void removeFvBoundary();

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

        //- Return cell face motion fluxes
        surfaceScalarField& phi()
        {
            return fvMesh::setPhi();
        }

        //- Return old-time cell volumes
        scalarField& V0()
        {
            return fvMesh::setV0();
        }

        //- Change topology of the mesh given the topological changes,
        //  change geometry and map the fields.
        void updateTopology(const polyTopoChange& changer);

        //- Reset time of last morph
        bool setMorphTimeIndex(const label newTimeIndex) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
