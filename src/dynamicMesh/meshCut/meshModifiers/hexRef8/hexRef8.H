/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    hexRef8

Description
    Splitting hexes into 8. Will abort if no hexes.

    Parallel: does communication to make sure faces across processor
    boundaries refine consistently. Will only force face to be refined,
    near-procpatch-cell does not get refined.

SourceFiles
    hexRef8.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef8_H
#define hexRef8_H

#include "labelList.H"
#include "face.H"
#include "className.H"
#include "labelHashSet.H"
#include "DynamicList.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class polyPatch;
class polyTopoChange;

/*---------------------------------------------------------------------------*\
                           Class hexRef8 Declaration
\*---------------------------------------------------------------------------*/

class hexRef8
{
    // Private data

        //- Reference to underlying mesh.
        const polyMesh& mesh_;

        //- Per cell the added cells (+ itself). Empty list if cell not
        //  refined.
        labelListList addedCells_;


    // Private Member Functions

        //- Synchronize cut faces across coupled patches
        void syncCoupledCutFaces(labelHashSet& cutFaces) const;

        //- Collect all cut boundaryEdges of pp. Returns these in patchFace
        //  and starting index in face (edge is f[index] to f[index+1])
        void collectCutEdges
        (
            const labelHashSet& cutEdges,
            const polyPatch& pp,
            DynamicList<label>& cutFaces,
            DynamicList<label>& cutIndex
        ) const;

        //- Synchronize cut edges across coupled patches. Returns true if
        //  anything changed (locally)
        bool syncCoupledCutEdges(labelHashSet& cutEdges) const;

        //- Calculate point-cell addressing. Master point first in list.
        void calcSortedCellPoints(labelListList& cellPoints) const;

        //- Get the two faces connected to edgeI on cellI in consistent order
        //  with e[0]
        void getEdgeFaces
        (
            const label cellI,
            const label edgeI,
            label&,
            label&
        ) const;

        //- Get the two edges using pointI on faceI. Such that e0 uses pointI
        //  and nextPointI
        void getFaceEdges
        (
            const label faceI,
            const label pointI,
            const label nextPointI,
            label& e0,
            label& e1
        ) const;

        //- Get patch and zone info
        void getFaceInfo
        (
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        ) const;

        //- Adds a face on top of existing faceI. Reverses if nessecary.
        void addFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Adds internal face from split of edge. No checks on reversal.
        void addInternalFace
        (
            polyTopoChange& meshMod,
            const label masterEdgeI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Modifies existing faceI for either new owner/neighbour or new face
        //  points. Reverses if nessecary.
        void modFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Look up new cell at cellI, pointI
        label newCell
        (
            const labelListList& cellPoints,
            const label cellI,
            const label pointI
        ) const;

        //- Get new owner and neighbour (in unspecified order) of pointI
        //  on faceI.
        void getFaceNeighbours
        (
            const labelListList& cellPoints,
            const label faceI,
            const label pointI,

            label& own,
            label& nei
        ) const;


        //- Disallow default bitwise copy construct
        hexRef8(const hexRef8&);

        //- Disallow default bitwise assignment
        void operator=(const hexRef8&);


public:

    //- Runtime type information
    ClassName("hexRef8");


    // Constructors

        //- Construct from mesh.
        hexRef8(const polyMesh& mesh);


    // Member Functions

        // Edit

            //- Insert refinement
            void setRefinement(const labelList& cells, polyTopoChange& meshMod);


        // Access

            //- From original to added cells. (includes original cell!)
            const labelListList& addedCells() const
            {
                return addedCells_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
