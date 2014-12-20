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
    mapAddedPolyMesh

Description
    Class containing mesh-to-mesh mapping information after a mesh addition
    where we add a mesh ('added mesh') to an old mesh, creating a new mesh.

    We choose to store mapping from the new mesh to the old and added meshes
    and if nessecary construct the inverse mapping. (not for the patches
    since patches might have been deleted). ?Might be a better idea to
    reverse this: store from old to current mesh and calculate if nessecary
    combined meshes?

    These maps are stored as a single labelList where a negative index is
    an index into the old mesh and a positive one into the added mesh. 
    (similar to face turning index in decomposePar)

    E.g.

    cellMap[i] < 0:

        oldCellI = -cellMap[I]-1;

    cellMap[i] > 0:

        addedCellI = cellMap[I]-1;

    cellMap[i] == 0:

        Illegal


    Note: Might need some more access functions or maybe some zone maps?


SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef mapAddedPolyMesh_H
#define mapAddedPolyMesh_H

//#include "refCount.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapAddedPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapAddedPolyMesh
{
    // Private data

        //- Old mesh points/face/cells
        label nOldPoints_;
        label nOldFaces_;
        label nOldCells_;

        //- Added mesh points/faces/cells
        label nAddedPoints_;
        label nAddedFaces_;
        label nAddedCells_;


        //- Old point map.
        //  Contains the point label on either old mesh or added mesh for all
        //  new points. Size of the list equals the size of new points
        labelList pointMap_;

        //- Old face map.
        //  Contains the face label on either old mesh or added mesh for all
        //  new faces. Size of the list equals the size of new faces
        labelList faceMap_;

        //- Old cell map.
        //  Contains the cell label on either old mesh or added mesh for all
        //  new cells. Size of the list equals the size of new cells
        labelList cellMap_;

        //- original mesh to new mesh patch map. -1 for deleted patches.
        labelList oldPatchMap_;

        //- added mesh to new mesh patch map. -1 for deleted patches.
        labelList addedPatchMap_;


    // Private member functions

        static labelList calcReverseMap
        (
            const labelList&,
            const bool getOldMesh,
            const label nOldElems
        );

public:

    // Constructors

        //- Construct null
        mapAddedPolyMesh();

        //- Construct from components
        mapAddedPolyMesh
        (
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            const label nAddedPoints,
            const label nAddedFaces,
            const label nAddedCells,
            const labelList& pointMap,
            const labelList& faceMap,
            const labelList& cellMap,
            const labelList& oldPatchMap,
            const labelList& addedPatchMap
        );

        //- Construct as copy
        mapAddedPolyMesh(const mapAddedPolyMesh&);


        // Helper functions for constructing.

            //- combine and invert map from old to new and from added to new
            static labelList combine
            (
                const label nNew,
                const labelList&,
                const labelList&
            );

            //- combine and invert map from old to new and from added to new.
            //  Difference with above one is that map from old to new is
            //  an identimap.
            static labelList combine
            (
                const label nNew,
                const label nOld,
                const labelList&
            );


    // Member Functions

        // Access

            // Old mesh data

                label nOldPoints() const
                {
                    return nOldPoints_;
                }

                label nOldFaces() const
                {
                    return nOldFaces_;
                }

                label nOldCells() const
                {
                    return nOldCells_;
                }

                //- From old patch index to new patch index or -1 if patch
                //  not present (since 0 size)
                const labelList& oldPatchMap() const
                {
                    return oldPatchMap_;
                }


            // Added mesh data

                label nAddedPoints() const
                {
                    return nAddedPoints_;
                }

                label nAddedFaces() const
                {
                    return nAddedFaces_;
                }

                label nAddedCells() const
                {
                    return nAddedCells_;
                }

                //- From added mesh patch index to new patch index or -1 if
                //  patch not present (since 0 size)
                const labelList& addedPatchMap() const
                {
                    return addedPatchMap_;
                }


            // Merged mesh data

                const labelList& pointMap() const
                {
                    return pointMap_;
                }

                const labelList& faceMap() const
                {
                    return faceMap_;
                }

                const labelList& cellMap() const
                {
                    return cellMap_;
                }


            //- From old mesh cells to new mesh cells.
            labelList oldCellMap() const;

            //- From added mesh cells to new mesh cells. -1 if cell not added.
            labelList addedCellMap() const;

            //- From old mesh faces to new mesh faces.
            labelList oldFaceMap() const;

            //- From added mesh faces to new mesh facess. -1 if face not added.
            labelList addedFaceMap() const;

            //- From old mesh points to new mesh points.
            labelList oldPointMap() const;

            //- From added mesh points to new mesh points. -1 if pt not added.
            labelList addedPointMap() const;


        // Edit

            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapAddedPolyMesh::updateMesh(const mapPolyMesh&)"
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
