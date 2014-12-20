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
    meshTools

Description
    Collection of static functions to do various simple mesh related things.

SourceFiles
    meshTools.C

\*---------------------------------------------------------------------------*/

#ifndef meshTools_H
#define meshTools_H

#include "label.H"
#include "vector.H"
#include "labelList.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "labelHashSet.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;


/*---------------------------------------------------------------------------*\
                           Class meshTools Declaration
\*---------------------------------------------------------------------------*/

class meshTools
{
    //- Bit identifiers for octants (p=plus, m=min e.g. plusXminYminZ)
    static const label mXmYmZ = 0;
    static const label pXmYmZ = 1;
    static const label mXpYmZ = 2;
    static const label pXpYmZ = 3;

    static const label mXmYpZ = 4;
    static const label pXmYpZ = 5;
    static const label mXpYpZ = 6;
    static const label pXpYpZ = 7;

    static const label mXmYmZMask = 1 << mXmYmZ;
    static const label pXmYmZMask = 1 << pXmYmZ;
    static const label mXpYmZMask = 1 << mXpYmZ;
    static const label pXpYmZMask = 1 << pXpYmZ;

    static const label mXmYpZMask = 1 << mXmYpZ;
    static const label pXmYpZMask = 1 << pXmYpZ;
    static const label mXpYpZMask = 1 << mXpYpZ;
    static const label pXpYpZMask = 1 << pXpYpZ;


public:

    // Static Functions

        //- Check if n is in same direction as normals of all faceLabels
        static bool visNormal
        (
            const vector& n,
            const vectorField& faceNormals,
            const labelList& faceLabels
        );

        //- Calculate point normals on a 'box' mesh (all edges aligned with
        //  coordinate axes)
        static vectorField calcBoxPointNormals(const primitivePatch& pp);

       //- Normalized edge vector
        static vector normEdgeVec(const primitiveMesh&, const label edgeI);

        //- Write obj representation of point
        static void writeOBJ
        (
            Ostream& os,
            const point& pt
        );


        //- Write obj representation of faces subset
        static void writeOBJ
        (
            Ostream& os,
            const faceList&,
            const pointField&,
            const labelList& faceLabels
        );

        //- Write obj representation of faces
        static void writeOBJ
        (
            Ostream& os,
            const faceList&,
            const pointField&
        );

        //- Write obj representation of cell subset
        static void writeOBJ
        (
            Ostream& os,
            const cellList&,
            const faceList&,
            const pointField&,
            const labelList& cellLabels
        );

        //- Is edge used by cell
        static bool edgeOnCell
        (
            const primitiveMesh&,
            const label cellI,
            const label edgeI
        );

        //- Is edge used by face
        static bool edgeOnFace
        (
            const primitiveMesh&,
            const label faceI,
            const label edgeI
        );

        //- Is face used by cell
        static bool faceOnCell
        (
            const primitiveMesh&,
            const label cellI,
            const label faceI
        );

        //- Return edge among candidates that uses the two vertices.
        static label findEdge
        (
            const edgeList& edges,
            const labelList& candidates,
            const label v0,
            const label v1
        );

        //- Return edge between two vertices. Returns -1 if no edge.
        static label findEdge
        (
            const primitiveMesh&,
            const label v0,
            const label v1
        );

        //- Return edge shared by two faces. Throws error if no edge found.
        static label getSharedEdge
        (
            const primitiveMesh&,
            const label f0,
            const label f1
        );

        //- Return face shared by two cells. Throws error if none found.
        static label getSharedFace
        (
            const primitiveMesh&,
            const label cell0,
            const label cell1
        );

        //- Get faces on cell using edgeI. Throws error if no two found.
        static void getEdgeFaces
        (
            const primitiveMesh&,
            const label cellI,
            const label edgeI,
            label& face0,
            label& face1
        );

        //- Return label of other edge (out of candidates edgeLabels)
        //  connected to vertex but not edgeI. Throws error if none found.
        static label otherEdge
        (
            const primitiveMesh&,
            const labelList& edgeLabels,
            const label edgeI,
            const label vertI
        );

        //- Return face on cell using edgeI but not faceI. Throws error 
        //  if none found.
        static label otherFace
        (
            const primitiveMesh&,
            const label cellI,
            const label faceI,
            const label edgeI
        );

        //- Return cell on other side of face. Throws error 
        //  if face not internal.
        static label otherCell
        (
            const primitiveMesh&,
            const label cellI,
            const label faceI
        );

        //- Returns label of edge nEdges away from startEdge (in the direction
        // of startVertI)
        static label walkFace
        (
            const primitiveMesh&,
            const label faceI,
            const label startEdgeI,
            const label startVertI,
            const label nEdges
        );


        //
        // Hex only functionality.
        //

        //- Given edge on hex find other 'parallel', non-connected edges.
        static void getParallelEdges
        (
            const primitiveMesh&,
            const label cellI,
            const label e0,
            label&,
            label&,
            label&
        );

        //- Given edge on hex find all 'parallel' (i.e. non-connected)
        //  edges and average direction of them
        static vector edgeToCutDir
        (
            const primitiveMesh&,
            const label cellI,
            const label edgeI
        );

        //- Reverse of edgeToCutDir: given direction find edge bundle and
        //  return one of them.
        static label cutDirToEdge
        (
            const primitiveMesh&,
            const label cellI,
            const vector& cutDir
        );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
