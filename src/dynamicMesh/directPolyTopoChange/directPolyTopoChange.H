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
    directPolyTopoChange

Description
    Simple mesh changes based on polyTopoChange syntax.

    Instead of recording changes and executing them all in one go (as does
    polyTopoChange) this class actually holds the current points/faces/cells
    and does the change immediately.
    It can be asked to compress out all unused points/faces/cells and
    renumber everything to be consistent.

    Note:
    - adding a face using non-existing cells causes all intermediate cells
    to be added. So always first add cells and then faces.
    (or set strict checking)
    - strict checking: -when adding faces owner and neighbour have to exist
        -no item can be removed more than once.
    - removed cell: cell set to 0 faces.
    - removed face: face set to 0 vertices.
    - removed point: coordinate set to greatPoint (GREAT,GREAT,GREAT).
    Note that this might give problems if this value is used already.
    To see if point is equal to above value we don't use == (which might give
    problems with roundoff error) but instead compare the individual component
    with >.
    - coupled patches: the reorderCoupledFaces routine (borrowed from
    the couplePatches utility) should reorder faces inside of the patches.
    It will be called automatically on the mesh(const polyBoundaryMesh&, ...)
    call but has to be called specifically if you are constructing your own
    patches.
    - maps (pointMap_, faceMap_, cellMap_) are from original mesh
    (or rather the order in which the items were added) to the new mesh.
    They are -1 for removed items. 

SourceFiles
    directPolyTopoChange.C
    directPolyTopoChangeCouples.C
    directPolyTopoChangeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef directPolyTopoChange_H
#define directPolyTopoChange_H

#include "autoPtr.H"
#include "DynamicList.H"
#include "labelList.H"
#include "IOobject.H"
#include "point.H"
#include "typeInfo.H"
#include "pointField.H"
#include "scalarField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class face;
class cell;
class polyMesh;
class polyMesh;
class Time;
class fileName;
class polyBoundaryMesh;
class polyPatch;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class directPolyTopoChange Declaration
\*---------------------------------------------------------------------------*/

class directPolyTopoChange
{
    // Private data

        //- Current point set
        DynamicList<point> points_;

        //- Current faceList
        DynamicList<face> faces_;

        //- Patch for every external face (-1 for internal faces)
        DynamicList<label> regions_;
        //- Owner for all faces
        DynamicList<label> faceOwner_;
        //- Neighbour for internal faces (-1 for external faces)
        DynamicList<label> faceNeighbour_;

        //- Current cell list
        DynamicList<cell> cells_;

        //- From original point order to compacted order
        DynamicList<label> pointMap_;
        DynamicList<label> faceMap_;
        DynamicList<label> cellMap_;

        //- Whether to allow referencing illegal points/cells/faces
        //  when adding/removing data.
        const bool strict_;


    // Private Member Functions

        //- Reorder contents of list according to map
        template<class T>
        static void reorder(const labelList& map, DynamicList<T>&);

        //- Reorder contents of list of list according to map
        template<class T>
        static void reorder(const labelList& map, List<DynamicList<T> >&);

        //- Renumber elements of list according to map
        static void renumber(const labelList&, DynamicList<label>&);

        //- Renumber & compact elements of list according to map
        static void renumberCompact(const labelList&, labelList&);

        //- Remove selected element from labelList
        static void remove(const label elemToRemove, labelList&);

        //- Is point removed?
        static bool pointRemoved(const point&);

        //- Do upper-triangular ordering and patch ordering.
        labelList getFaceOrder(const label nPatches) const;

        //- Compact and reorder faces according to map. Optionally check if
        //  compacting faces creates illegal cells.
        void reorderCompactFaces
        (
            const label newSize,
            const labelList& oldToNew,
            const bool checkCells
        );

        //- Sizes of patches and number of internal faces. patchSizes/starts
        //  are either number of old patches or the max. region number + 1,
        //  whichever is bigger.
        void calcPatchSizes
        (
            const label nOldPatches,
            label& nInternalFaces,
            labelList& patchSizes,
            labelList& patchStarts
        ) const;

        //- Create list of patches from patch dictionaries
        static List<polyPatch*> createPatches
        (
            const wordList& patchNames,
            const PtrList<dictionary>& patchDicts,
            const label nInternalFaces,
            const labelList& patchSizes,
            const polyBoundaryMesh& boundaryMesh
        );


        // Coupling

            //- Get face centre for slice of faces.
            pointField calcFaceCentres(const label start, const label size)
             const;

            //- Get coordinate of f[0] for slice of faces
            pointField getAnchorPoints(const label start, const label size)
             const;

            // Get local typical dimension and tolerance from that. Currently
            // max of distance from centre to any of the face points.
            scalarField calcFaceTol
            (
                const label start,
                const label size,
                const pointField& faceCentres
            ) const;

            //- Calculate anticlockwise rotation for f[0] to align with anchor
            label getRotation
            (
                const face& f,
                const point& anchor,
                const scalar tol
            ) const;

            //- Rotate face by number of positions
            static face rotateFace(const face& f, const label nPos);


            //- Send ordering information on processorPolyPatch
            void initOrderProcessorPolyPatch
            (
                const label neighbourProcNo,
                const label start,
                const label size
            ) const;

            //- Send ordering information on cyclicPolyPatch
            void initOrderCyclicPolyPatch
            (
                const label start,
                const label size
            ) const;

            //- Receive and calulate ordering
            bool orderProcessorPolyPatch
            (
                const label neighbourProcNo,
                const label start,
                const label size,

                labelList& faceMap,
                labelList& rotation
            ) const;

            //- Receive and calulate ordering
            bool orderCyclicPolyPatch
            (
                const label start,
                const label size,

                labelList& faceMap,
                labelList& rotation
            ) const;

            //- Do all coupled patch face reordering
            void reorderCoupledFaces
            (
                const PtrList<dictionary>& oldPatchDicts,
                const labelList& patchStarts,
                const labelList& patchSizes
            );

public:

    //- Runtime type information
    ClassName("directPolyTopoChange");


    // Static data members

        //- Value of deleted point
        static const point greatPoint;


    // Constructors

        //- Construct null.
        directPolyTopoChange(const bool strict = true);

        //- Construct from mesh. Adds all points/face/cells from mesh.
        directPolyTopoChange(const polyMesh& mesh, const bool strict = true);


    // Member Functions

        // Access

            //- Points. Shrunk after constructing mesh (or calling of compact())
            const DynamicList<point>& points() const
            {
                return points_;
            }

            const DynamicList<face>& faces() const
            {
                return faces_;
            }

            const DynamicList<label>& regions()const
            {
                return regions_;
            }

            const DynamicList<label>& faceOwner()const
            {
                return faceOwner_;
            }

            const DynamicList<label>& faceNeighbour()const
            {
                return faceNeighbour_;
            }

            const DynamicList<cell>& cells() const
            {
                return cells_;
            }

            const DynamicList<label>& pointMap()const
            {
                return pointMap_;
            }

            const DynamicList<label>& faceMap()const
            {
                return faceMap_;
            }

            const DynamicList<label>& cellMap()const
            {
                return cellMap_;
            }


        // Edit

            //- Clear all storage
            void clear();

            //- Add all points/faces/cells of mesh. Additional offset for patch
            //  ids.
            void addMesh(const polyMesh&, const label patchOffset = 0);


            //- Add point. Return new point label.
            label addPoint(const point&);

            //- Modify coordinate.
            void modifyPoint(const label, const point&);

            //- Bulk modify coordinates
            void movePoints(const pointField&);

            //- Remove point.
            void removePoint(const label);


            //- Add face to cells. Return new face label.
            label addFace
            (
                const face&,
                const label own,
                const label nei,
                const label patchI
            );

            //- Modify vertices or cell of face.
            void modifyFace
            (
                const face&,
                const label faceI,
                const label own,
                const label nei,
                const label patchI
            );

            //- Remove face.
            void removeFace(const label);


            //- Add cell. Return new cell label.
            label addCell();

            //- Remove cell.
            void removeCell(const label);


            //- Remove all unused/removed points/faces/cells and update
            //  face ordering.
            void compact();


        // Other

            //- Construct mesh without patches. Use addPatches to add patches.
            //  Clears all data (faces_, cells_, etc) but the maps and regions_.
            autoPtr<polyMesh> mesh(const IOobject& io);

            //- Helper function to get patch data in form of dictionaries
            static PtrList<dictionary> getPatchDicts(const polyBoundaryMesh&);

            //- Construct mesh with patch names/types from old mesh.
            //  Adds old patches (and any new ones) and does
            //  reorderCoupledFaces.
            //  Clears all data (faces_, cells_, etc) but the maps.
            autoPtr<polyMesh> mesh
            (
                const wordList& oldPatchNames,
                const PtrList<dictionary>& oldPatchDicts,
                const IOobject& io
            );

            //- Inplace changes mesh with patch names/types from old mesh.
            //  Adds old patches (and any new ones) and does
            //  reorderCoupledFaces.
            //  Clears all data (faces_, cells_, etc) but the maps.
            void makeMesh
            (
                const wordList& oldPatchNames,
                const PtrList<dictionary>& oldPatchDicts,
                polyMesh& mesh
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "directPolyTopoChangeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
