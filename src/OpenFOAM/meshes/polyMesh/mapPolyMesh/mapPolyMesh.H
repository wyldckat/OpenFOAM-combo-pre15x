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
    mapPolyMesh

Description
    Class containing mesh-to-mesh mapping information after a change
    in polyMesh topology.

SourceFiles
    mapPolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mapPolyMesh_H
#define mapPolyMesh_H

#include "labelList.H"
#include "objectMap.H"
#include "pointField.H"
#include "labelHashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapPolyMesh
:
    public refCount
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Number of old live points
        label nOldPoints_;

        //- Number of old live faces
        label nOldFaces_;

        //- Number of old live cells
        label nOldCells_;

        //- Old point map.
        //  Contains the old point label for all new points.
        //  - for preserved points this is the old point label.
        //  - for added points this is the master point ID
        //  - for points added with no master, this is -1
        //  Size of the list equals the size of new points
        labelList pointMap_;

        //- Old face map.
        //  Contains a list of old face labels for every new face.
        //  Size of the list equals the number of new faces
        //  - for preserved faces this is the old face label.
        //  - for faces added from faces this is the master face ID
        //  - for faces added with no master, this is -1
        //  - for faces added from points or edges, this is -1
        labelList faceMap_;

        //- Faces inflated from points
        List<objectMap> facesFromPointsMap_;

        //- Faces inflated from edges
        List<objectMap> facesFromEdgesMap_;

        //- Old cell map.
        //  Contains old cell label for all preserved cells.
        //  Size of the list equals the number or preserved cells
        labelList cellMap_;

        //- Cells inflated from points
        List<objectMap> cellsFromPointsMap_;

        //- Cells inflated from edges
        List<objectMap> cellsFromEdgesMap_;

        //- Cells inflated from faces
        List<objectMap> cellsFromFacesMap_;

        //- Reverse point map
        labelList reversePointMap_;

        //- Reverse face map
        labelList reverseFaceMap_;

        //- Reverse cell map
        labelList reverseCellMap_;

        //- Map of flipped face flux faces
        labelHashSet flipFaceFlux_;

        //- Patch mesh point renumbering
        labelListList patchPointMap_;

        //- Point zone renumbering
        //  For every preserved point in zone give the old position.
        //  For added points, the index is set to -1
        labelListList pointZoneMap_;

        //- Face zone point renumbering
        //  For every preserved point in zone give the old position.
        //  For added points, the index is set to -1
        labelListList faceZonePointMap_;

        //- Face zone face renumbering
        //  For every preserved face in zone give the old position.
        //  For added faces, the index is set to -1
        labelListList faceZoneFaceMap_;

        //- Cell zone renumbering
        //  For every preserved cell in zone give the old position.
        //  For added cells, the index is set to -1
        labelListList cellZoneMap_;

        //- Pre-motion point positions.
        //  This specifies the correct way of blowing up zero-volume objects
        pointField preMotionPoints_;

        //- List of the old patch sizes
        labelList oldPatchSizes_;

        //- List of the old patch start labels
        labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        labelList oldPatchNMeshPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mapPolyMesh(const mapPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const mapPolyMesh&);

        // Private member functions to calculate demand driven data


public:

    // Constructors

        //- Construct from components
        mapPolyMesh
        (
            const polyMesh& mesh,
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            const labelList& pointMap,
            const labelList& faceMap,
            const List<objectMap>& facesFromPoints,
            const List<objectMap>& facesFromEdges,
            const labelList& cellMap,
            const List<objectMap>& cellsFromPoints,
            const List<objectMap>& cellsFromEdges,
            const List<objectMap>& cellsFromFaces,
            const labelList& reversePointMap,
            const labelList& reverseFaceMap,
            const labelList& reverseCellMap,
            const labelHashSet& flipFaceFlux,
            const labelListList& patchPointMap,
            const labelListList& pointZoneMap,
            const labelListList& faceZonePointMap,
            const labelListList& faceZoneFaceMap,
            const labelListList& cellZoneMap,
            const pointField& preMotionPoints,
            const labelList& oldPatchStarts,
            const labelList& oldPatchNMeshPoints
        );


    // Member Functions

        // Access

            //- Return polyMesh
            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Number of old points
            label nOldPoints() const
            {
                return nOldPoints_;
            }

            //- Number of old internal faces
            label nOldInternalFaces() const
            {
                return oldPatchStarts_[0];
            }

            //- Number of old faces
            label nOldFaces() const
            {
                return nOldFaces_;
            }

            //- Number of old cells
            label nOldCells() const
            {
                return nOldCells_;
            }

            //- Old point map.
            //  Contains the old point label for all new points.
            //  For preserved points this is the old point label.
            //  For added points this is the master point ID
            const labelList& pointMap() const
            {
                return pointMap_;
            }

            //- Old face map.
            //  Contains a list of old face labels for every new face.
            //  Warning: this map contains invalid entries for new faces
            const labelList& faceMap() const
            {
                return faceMap_;
            }

            //- Faces inflated from points
            const List<objectMap>& facesFromPointsMap() const
            {
                return facesFromPointsMap_;
            }

            //- Faces inflated from edges
            const List<objectMap>& facesFromEdgesMap() const
            {
                return facesFromEdgesMap_;
            }

            //- Old cell map.
            //  Contains old cell label for all preserved cells.
            const labelList& cellMap() const
            {
                return cellMap_;
            }

            //- Cells inflated from points
            const List<objectMap>& cellsFromPointsMap() const
            {
                return cellsFromPointsMap_;
            }

            //- Cells inflated from edges
            const List<objectMap>& cellsFromEdgesMap() const
            {
                return cellsFromEdgesMap_;
            }

            //- Cells inflated from faces
            const List<objectMap>& cellsFromFacesMap() const
            {
                return cellsFromFacesMap_;
            }


            // Reverse maps

                //- Reverse point map
                //  Contains new point label for all old and added points
                const labelList& reversePointMap() const
                {
                    return reversePointMap_;
                }


                //- Reverse face map
                //  Contains new face label for all old and added faces
                const labelList& reverseFaceMap() const
                {
                    return reverseFaceMap_;
                }

                //- Reverse cell map
                //  Contains new cell label for all old and added cells
                const labelList& reverseCellMap() const
                {
                    return reverseCellMap_;
                }


            //- Map of flipped face flux faces
            const labelHashSet& flipFaceFlux() const
            {
                return flipFaceFlux_;
            }

                //- Patch point renumbering
                //  For every preserved point on a patch give the old position.
                //  For added points, the index is set to -1
                const labelListList& patchPointMap() const
                {
                    return patchPointMap_;
                }


            // Zone mapping

                //- Point zone renumbering
                //  For every preserved point in zone give the old position.
                //  For added points, the index is set to -1
                const labelListList& pointZoneMap() const
                {
                    return pointZoneMap_;
                }

                //- Face zone point renumbering
                //  For every preserved point in zone give the old position.
                //  For added points, the index is set to -1
                const labelListList& faceZonePointMap() const
                {
                    return faceZonePointMap_;
                }

                //- Face zone face renumbering
                //  For every preserved face in zone give the old position.
                //  For added faces, the index is set to -1
                const labelListList& faceZoneFaceMap() const
                {
                    return faceZoneFaceMap_;
                }

                //- Cell zone renumbering
                //  For every preserved cell in zone give the old position.
                //  For added cells, the index is set to -1
                const labelListList& cellZoneMap() const
                {
                    return cellZoneMap_;
                }

            //- Pre-motion point positions.
            //  This specifies the correct way of blowing up
            //  zero-volume objects
            const pointField& preMotionPoints() const 
            {
                return preMotionPoints_;
            }

            //- Return list of the old patch sizes
            const labelList& oldPatchSizes() const
            {
                return oldPatchSizes_;
            }

            //- Return list of the old patch start labels
            const labelList& oldPatchStarts() const
            {
                return oldPatchStarts_;
            }

            //- Return numbers of mesh points per old patch
            const labelList& oldPatchNMeshPoints() const
            {
                return oldPatchNMeshPoints_;
            }

            //- Modify for reordered faces. Is very limited 'updateTopology';
            //  only allows for reshuffling of faces inside of a patch.
            //  Gets passed faceMap (old to new).
            void reorderPatchFaces(const labelList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
