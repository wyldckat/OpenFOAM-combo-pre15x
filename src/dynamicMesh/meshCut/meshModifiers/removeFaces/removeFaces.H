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
    removeFaces

Description
    Given list of faces to remove insert all the topology changes.

SourceFiles
    removeFaces.C

\*---------------------------------------------------------------------------*/

#ifndef removeFaces_H
#define removeFaces_H

#include "physicalConstants.H"
#include "labelList.H"
#include "Map.H"
#include "labelHashSet.H"
#include "boolList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class primitiveMesh;
class polyMesh;
class polyTopoChange;
class edge;
class face;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class removeFaces Declaration
\*---------------------------------------------------------------------------*/

class removeFaces
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Minimum cos of angle between faces to be considered aligned.
        //  (and thus candidates for merging)
        const scalar minCos_;


    // Private Static Functions

        //- Next index in list
        static label nexti(const label, const label index);

        //- Previous index in list
        static label previ(const label, const label index);

        //- Return true if f0 and f1 (with shared edge e) use the points
        //  of the edge in reverse order.
        static bool sameFaceOrdering
        (
            const edge& e,
            const face& f0,
            const label f0Start,    // index of edge start on f0
            const face& f1,
            const label f1Start     // index of edge start on f1
        );


    // Private Member Functions

        // Detection of valid coarsening patterns

            //- Get cos of angle between two faces shared by edge.
            //  If checkOrder it checks the order of shared face points and
            //  reverses the angle if nessecary.
            scalar faceCos
            (
                const label face0I,
                const label face1I,
                const bool checkOrder,
                const label sharedEdgeI
            ) const;

        // Detection of coarsening influence

            //- Get patch, zone info for faceI
            void getFaceInfo
            (
                const label faceI,

                label& patchID,
                label& zoneID,
                label& zoneFlip
            ) const;

            //- Count number of remaining faces using edgeI. Returns true if
            //  only two remaining and sets face0I, face1I to those two.
            bool getUnremovedFaces
            (
                const labelHashSet& facesToRemove,
                const label edgeI,

                label& face0I,
                label& face1I
            ) const;

            //- Count number of remaining edges using vertI. Returns true if 
            //  only two remaining and sets edge0I, edge1I to those two.
            bool getUnremovedEdges
            (
                const Map<edge>& edgesToRemove,
                const label vertI,

                label& edge0I,
                label& edge1I
            ) const;

            //- Checks vertI for if it can be removed. (so if only two edges
            //  using it and these two are in line).
            bool canBeRemoved
            (
                const Map<edge>& edgesToRemove,
                const label vertI
            ) const;

            //- Get all affected faces (including faces marked for removal)
            labelHashSet getFacesAffected
            (
                const labelList& cellRegion,
                const labelList& cellRegionMaster,
                const labelHashSet& facesToRemove,
                const Map<edge>& edgesToRemove,
                const labelHashSet& pointsToRemove
            ) const;


        // Topological changes

            //- Return face with all pointsToRemove removed.
            face filterFace(const labelHashSet& pointsToRemove, const label)
                const;

            //- Combines two faces, leaving out points on the edge they share
            //  and all points in pointsToRemove. Ordering of points is as
            //  face0.
            face mergeFace
            (
                const labelHashSet& pointsToRemove,
                const label face0I,
                const label face1I,
                const label sharedEdgeI
            ) const;

            //- Wrapper for meshMod.modifyFace. Reverses face if own>nei.
            void modFace
            (
                const face& f,
                const label masterFaceID,
                const label own,
                const label nei,
                const bool flipFaceFlux,
                const label newPatchID,
                const bool removeFromZone,
                const label zoneID,
                const bool zoneFlip,

                polyTopoChange& meshMod
            ) const;


        //- Change elements in cellRegion that are oldRegion to newRegion.
        //  Recurses to cell neighbours.
        void changeCellRegion
        (
            const label cellI,
            const label oldRegion,
            const label newRegion,
            labelList& cellRegion
        ) const;



        //- Disallow default bitwise copy construct
        removeFaces(const removeFaces&);

        //- Disallow default bitwise assignment
        void operator=(const removeFaces&);


public:

    //- Runtime type information
    ClassName("removeFaces");


    // Constructors

        //- Construct from mesh and min cos of angle for faces/edges to be
        //  considered aligned.
        removeFaces
        (
            const polyMesh& mesh,
            const scalar minCos = Foam::cos(10.0/180.0 * physicalConstant::pi)
        );


    // Member Functions

        // Access

            scalar minCos() const
            {
                return minCos_;
            }

        // Helper

            //- Given set of faces to pierce calculates:
            //  - region for connected cells
            //  - mastercell for each region
            //  - new set of faces which contains input set + additional ones
            //    where cells on both sides would have same mastercell.
            //  Returns number of regions.
            label compatibleRemoves
            (
                const labelList& inPiercedFaces,
                labelList& cellRegion,
                labelList& cellRegionMaster,
                labelList& outPiercedFaces
            ) const;


        // Topology changes

            //- Play commands into polyTopoChange to remove faces.
            void setRefinement
            (
                const labelList& piercedFaces,
                const labelList& cellRegion,
                const labelList& cellRegionMaster,
                polyTopoChange&
            ) const;

            //- Force recalculation of locally stored data on topological change
            void updateTopology(const mapPolyMesh&)
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
