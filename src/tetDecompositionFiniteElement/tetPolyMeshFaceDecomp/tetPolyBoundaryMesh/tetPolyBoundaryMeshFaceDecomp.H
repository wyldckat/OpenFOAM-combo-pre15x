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
    tetPolyBoundaryMeshFaceDecomp

Description

SourceFiles
    tetPolyBoundaryMeshFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyBoundaryMeshFaceDecomp_H
#define tetPolyBoundaryMeshFaceDecomp_H

#include "tetPolyPatchFaceDecompList.H"
#include "faceListFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMeshFaceDecomp;
class polyBoundaryMesh;
class globalProcessorTetPolyPatchFaceDecomp;

/*---------------------------------------------------------------------------*\
                 Class tetPolyBoundaryMeshFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyBoundaryMeshFaceDecomp
:
    public tetPolyPatchFaceDecompList
{
    // Private data

        //- Reference to mesh
        const tetPolyMeshFaceDecomp& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyBoundaryMeshFaceDecomp(const tetPolyBoundaryMeshFaceDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyBoundaryMeshFaceDecomp&);


public:

    // Constructors

        //- Construct from polyBoundaryMesh
        tetPolyBoundaryMeshFaceDecomp
        (
            const tetPolyMeshFaceDecomp&,
            const polyBoundaryMesh&
        );


    // Member functions

        //- Return the mesh reference
        const tetPolyMeshFaceDecomp& mesh() const
        {
            return mesh_;
        }

        //- Return reference to shared processor point patch
        const tetPolyPatchFaceDecomp& globalPointPatch() const;

        //- Return boundary faces per patch.  All are triangles
        faceListList boundaryTriFaces() const;

        // Update topology
        void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
