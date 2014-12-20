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
    - coupled patches: the reorderCoupledPatches routine (borrowed from
    the couplePatches utility) should reorder faces inside of the patches.
    It will be called automatically on the mesh(const polyBoundaryMesh&, ...)
    call but has to be called specifically if you are constructing your own
    patches.
    - maps (pointMap_, faceMap_, cellMap_) are from original mesh
    (or rather the order in which the items were added) to the new mesh.
    They are -1 for removed items. 

SourceFiles
    directPolyTopoChange.C
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class face;
class cell;
class polyMesh;
class morphMesh;
class Time;
class fileName;
class polyBoundaryMesh;

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

        //- Sizes of patches and number of internal faces.
        void calcPatchSizes(label& nInternalFaces, labelList& patchSizes) const;


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

            //- Add all points/faces/cells of mesh.
            void addMesh(const polyMesh&);

            //- Add point. Return new point label.
            label addPoint(const point&);

            //- Modify coordinate.
            void modifyPoint(const label, const point&);

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
            //  (Call reorderCoupledPatches to update cyclics/procBoundaries)
            autoPtr<morphMesh> mesh(const IOobject& io);

            //- Sort faces on coupled patches to be consistent.
            //  Does parallel communication.
            void reorderCoupledPatches(morphMesh& mesh);

            //- Construct mesh with patch names/types from old mesh.
            //  Adds old patches (and any new ones) and does
            //  reorderCoupledPatches.
            //  Clears all data (faces_, cells_, etc) but the maps.
            autoPtr<morphMesh> mesh
            (
                const wordList& oldPatchNames,
                const wordList& oldPatchTypes,
                const IOobject& io
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
