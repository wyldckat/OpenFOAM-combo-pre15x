// The FOAM Project // File: meshTriangulation.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   meshTriangulation
   \\  /           | Family: triSurface
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2005 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    meshTriangulation

DESCRIPTION
    Triangulation (in compact point numbering) of selected mesh faces.
    Uses faceTriangulation so handles non-convex faces.

    Is the triangulated faces after construction. Note that this will usually
    be multiply connected.

C++
    listOfSourceFiles
        meshTriangulation.C
    endListOfSourceFiles

AUTHOR
    Mattijs Janssens.
*/
// ------------------------------------------------------------------------- //

#ifndef meshTriangulation_H
#define meshTriangulation_H

#include "triSurface.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class meshTriangulation Declaration
\*---------------------------------------------------------------------------*/

class meshTriangulation
:
    public triSurface
{
    // Private data

        //- Number of triangles in this that are internal to the surface.
        label nInternalFaces_;

        //- From triangle to mesh face
        labelList faceMap_;

    // Private Member Functions

        static bool isInternalFace
        (
            const primitiveMesh& mesh,
            const boolList& includedCell,
            const label faceI
        );


public:

    ClassName("meshTriangulation");


    // Constructors

        //- Construct null
        meshTriangulation();

        //- Construct from selected mesh cell and region number to be used
        //  for triangles resulting from internal faces. (all boundary triangles
        //  get polyMesh patch id)
        meshTriangulation
        (
            const polyMesh&,
            const label internalFacesPatch,
            const boolList& includedCell
        );


    // Member Functions

        //- number of triangles in *this which are internal to the surface
        label nInternalFaces() const
        {
            return nInternalFaces_;
        }

        //- from triangle to mesh face
        const labelList& faceMap() const
        {
            return faceMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
