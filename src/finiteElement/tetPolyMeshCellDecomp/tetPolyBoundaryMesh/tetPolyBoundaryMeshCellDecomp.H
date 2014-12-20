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
    tetPolyBoundaryMeshCellDecomp

Description

SourceFiles
    tetPolyBoundaryMeshCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyBoundaryMeshCellDecomp_H
#define tetPolyBoundaryMeshCellDecomp_H

#include "tetPolyPatchCellDecompList.H"
#include "faceListFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMeshCellDecomp;
class polyBoundaryMesh;
class globalProcessorTetPolyPatchCellDecomp;

/*---------------------------------------------------------------------------*\
                  Class tetPolyBoundaryMeshCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyBoundaryMeshCellDecomp
:
    public tetPolyPatchCellDecompList
{
    // Private data

        //- Reference to mesh
        const tetPolyMeshCellDecomp& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyBoundaryMeshCellDecomp(const tetPolyBoundaryMeshCellDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyBoundaryMeshCellDecomp&);


public:

    // Constructors

        //- Construct from polyBoundaryMesh
        tetPolyBoundaryMeshCellDecomp
        (
            const tetPolyMeshCellDecomp&,
            const polyBoundaryMesh&
        );


    // Member functions

        //- Return the mesh reference
        const tetPolyMeshCellDecomp& mesh() const
        {
            return mesh_;
        }

        //- Return reference to shared processor point patch
        const tetPolyPatchCellDecomp& globalPointPatch() const;

        //- Return boundary faces per patch.  All are triangles
        faceListList boundaryTriFaces() const;

        // Update topology
        void updateTopology();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
