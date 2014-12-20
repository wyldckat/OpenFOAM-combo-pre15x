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
    repatchPolyMesh

Description
    A mesh which allows changes in the patch distribution of the
    boundary faces.  The change in patching is set using changePatchID. For a
    boundary face, a new patch ID is given. 

SourceFiles
    repatchPolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef repatchPolyMesh_H
#define repatchPolyMesh_H

#include "polyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class repatchPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class repatchPolyMesh
:
    public polyMesh
{
    // Private data

        //- Topological change to accumulated all mesh changes
        polyTopoChange meshMod_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        repatchPolyMesh(const repatchPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const repatchPolyMesh&);


public:

    //- Runtime type information
    TypeName("repatchPolyMesh");


    // Constructors

        //- Construct from IOobject
        explicit repatchPolyMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addPatches() member function
        repatchPolyMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Construct from cell shapes
        repatchPolyMesh
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

        virtual ~repatchPolyMesh()
        {}


    // Member Functions

        //- Change patches.
        void changePatches(const List<polyPatch*>& patches)
        {
            polyMesh::removeBoundary();
            polyMesh::addPatches(patches);
        }

        //- Change patch ID for a boundary face. Note: patchID should be in new
        //  numbering.
        void changePatchID
        (
            const label faceID,
            const label patchID
        );

        //- Set zone ID for a face
        void setFaceZone
        (
            const label faceID,
            const label zoneID,
            const bool zoneFlip
        );

        //- Change anchor point (zero'th point of face) for a boundary face.
        void changeAnchorPoint
        (
            const label faceID,
            const label fp
        );


        //- Re-patch the mesh
        void repatch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
