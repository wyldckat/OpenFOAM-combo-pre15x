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
    attachPolyMesh

Description
    This class is derived from polyMesh and serves as a tool for
    statically connecting pieces of a mesh by executing the mesh
    modifiers and cleaning the mesh.

    The idea is that a mesh can be built from pieces and put together
    using various mesh modifiers (mainly sliding interfaces) which are
    not needed during the run.  Therefore, once the mesh is assembled
    and mesh modification triggered, the newly created point, face and
    cell zones can be cleared together with the mesh modifiers thus
    creating a singly connected static mesh.

    Note:
    All point, face and cell zoning will be lost!  Do it after
    attaching the parts of the mesh, as the point, face and cell
    numbering changes.

\*---------------------------------------------------------------------------*/

#ifndef attachPolyMesh_H
#define attachPolyMesh_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class attachPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class attachPolyMesh
:
    public polyMesh
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        attachPolyMesh(const attachPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const attachPolyMesh&);


public:

    // Constructors

        //- Construct from IOobject
        attachPolyMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addPatches() member function
        attachPolyMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Construct from cell shapes
        attachPolyMesh
        (
            const IOobject& io,
            const pointField& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const wordList& boundaryPatchTypes,
            const word& defaultBoundaryPatchType,
            const wordList& boundaryPatchPhysicalTypes
        );


    // Destructor

        virtual ~attachPolyMesh();


    // Member Functions

        //- Attach mesh. By default filter out empty patches.
        void attach(const bool removeEmptyPatches = true);


        //- Write after removing files from old mesh
        virtual bool write
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;


        bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
