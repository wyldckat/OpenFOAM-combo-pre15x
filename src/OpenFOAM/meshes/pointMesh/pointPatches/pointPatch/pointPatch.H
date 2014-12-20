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
    pointPatch

Description
    Basic pointPatch represents a set of points from the mesh.

SourceFiles
    pointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatch_H
#define pointPatch_H

#include "labelList.H"
#include "vectorField.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

class pointBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class pointPatch Declaration
\*---------------------------------------------------------------------------*/

class pointPatch
{
    // Private data

        //- Reference to boundary mesh
        const pointBoundaryMesh& boundaryMesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointPatch(const pointPatch&);

        //- Disallow default bitwise assignment
        void operator=(const pointPatch&);


public:

    //- Runtime type information
    TypeName("basePatch");


    // Constructor

        pointPatch
        (
            const pointBoundaryMesh& bm
        )
        :
            boundaryMesh_(bm)
        {}


    // Destructor

        virtual ~pointPatch()
        {}


    // Member Functions

        //- Return name
        virtual const word& name() const = 0;

        //- Return size
        virtual label size() const = 0;

        //- Return number of faces
        virtual label nFaces() const = 0;

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const = 0;

        //- Return boundaryMesh reference
        const pointBoundaryMesh& boundaryMesh() const
        {
            return boundaryMesh_;
        }

        //- Return true if this patch field is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const = 0;

        //- Return mesh points
        virtual const vectorField& localPoints() const = 0;

        //- Return  point normals
        virtual const vectorField& pointNormals() const = 0;

        //- Face decomposition into triangle
        //  Used for implementation of FEM boundary conditions
        virtual triFaceList faceTriangles(const label faceID) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
