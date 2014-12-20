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
    polyMesh

Description
    Moving poly mesh with topological changes.
    Supports general polyhedral cells.

SourceFiles
    polyMesh.C
    polyMeshIO.C
    calcPolyMeshFaceCells.C
    clearPolyMesh.C
    createPolyMesh.C
    polyMeshCorrectFaceDirections.C
    polyMeshMorph.C
    polyMeshPurgePolygons.C
    polyMeshUpdateTopology.C
    polyShapeMeshTools.C
    zipUpCells.C

\*---------------------------------------------------------------------------*/

#ifndef polyMesh_H
#define polyMesh_H

#include "objectRegistry.H"
#include "primitiveMesh.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "cellShapeList.H"
#include "pointIOField.H"
#include "faceIOList.H"
#include "cellIOList.H"
#include "polyBoundaryMesh.H"
#include "pointZoneMesh.H"
#include "faceZoneMesh.H"
#include "cellZoneMesh.H"
#include "polyMeshMorphEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class parallelInfo;

/*---------------------------------------------------------------------------*\
                          Class polyMesh Declaration
\*---------------------------------------------------------------------------*/

class polyMesh
:
    public objectRegistry,
    public primitiveMesh
{

public:

    // Public data types

        //- Enumeration defining the state of the mesh after a read update.
        //  Used for post-processing applications, where the mesh
        //  needs to update based on the files written in time
        //  directores
        enum readUpdateState
        {
            UNCHANGED,
            POINTS_MOVED,
            TOPO_CHANGE,
            TOPO_PATCH_CHANGE
        };


private:

    // Permanent data

        // Primitive mesh data

            //- Points
            pointIOField points_;

            //- Faces
            faceIOList faces_;

            //- Cells
            cellIOList cells_;

            //- Boundary mesh
            mutable polyBoundaryMesh boundary_;


            // Zoning information

                //- Point zones
                pointZoneMesh pointZones_;

                //- Face zones
                faceZoneMesh faceZones_;

                //- Cell zones
                cellZoneMesh cellZones_;


            //- Parallel info
            mutable parallelInfo* parallelDataPtr_;

            //- Topology morphing engine
            mutable polyMeshMorphEngine* morphEnginePtr_;


        // Minimal derived data

            //- Face owner
            mutable labelList* allOwnerPtr_;

            //- Face neighbour
            mutable labelList* allNeighbourPtr_;


        // Other mesh related data

            // Mesh motion

                //- Is the mesh moving
                bool moving_;

                //- Current time index for mesh motion
                mutable label curMotionTimeIndex_;

                //- Old points (for the last mesh motion)
                mutable pointField* oldPointsPtr_;

            // Topological change

                //- Is the mesh morphing topology
                mutable bool morphing_;

                //- Current time index for mesh motion
                mutable label curMorphTimeIndex_;

                //- Topology morphing map for current topological change
                mutable mapPolyMesh* morphMap_;


    // Private member functions

        //- Disallow construct as copy
        polyMesh(const polyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const polyMesh&);

        // Private member functions to calculate demand driven data

            //- Calculate cell-face addressing
            void calcFaceCells() const;

            //- Clear cell-face addressing
            void clearFaceCells() const;

            //- Calculate the cell shapes from the primitive
            //  polyhedral information
            void calcCellShapes() const;


        // Helper functions for constructor from cell shapes
    
            labelListList cellShapePointCells(const cellShapeList&) const;

            labelList facePatchFaceCells
            (
                const faceList& patchFaces,
                const labelListList& pointCells,
                const faceListList& cellsFaceShapes,
                const label patchID
            ) const;


        // Topology changes

            //- Change topology of the mesh given a topological change
            void morph(const polyTopoChange&);

            //- Rotate a face nShift positions in anticlockwise direction
            static face rotateFace(const face& f, const label nShift);

            //- Do everything nessecary to reorder faces on coupled patches
            //  consistently. Return true if anything changed. Updates
            //  all topological data and morphMap.
            bool reorderCoupledPatches
            (
                const polyTopoChange& ref,
                mapPolyMesh& morphMap
            );


protected:

        // Special mesh operations

            //- Update mesh topology using explicitly provided mesh changes.
            //  This is incompatible with morphEngine/meshModifiers mechanism
            //  and this function will FatalError if there is a morphEngine.
            //  This function is only to be used if you cannot do your changes
            //  through meshModifiers.
            void updateTopology(const polyTopoChange& ref);

            //- Reset morph time index
            bool setMorphTimeIndex(const label newTimeIndex) const;

            //- Remove boundary patches
            void removeBoundary();

            //- Remove point zones
            void removePointZones();

            //- Remove face zones
            void removeFaceZones();

            //- Remove cell zones
            void removeCellZones();

            //- Remove mesh modifiers
            void removeMeshModifiers();

            //- Set the instance for mesh files
            void setInstance(const fileName&);


public:

    //- Runtime type information
    TypeName("polyMesh");

    //- Additional topology morph debug switch
    static int morphDebug;

    //- Return the default region name
    static word defaultRegion;

    //- Return the mesh sub-directory name (usually "polyMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from IOobject
        explicit polyMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addPatches() member function
        polyMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Construct from cell shapes
        polyMesh
        (
            const IOobject& io,
            const pointField& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const wordList& boundaryPatchTypes,
            const word& defaultBoundaryPatchType,
            const wordList& boundaryPatchPhysicalTypes
        );


    // Destructor

        virtual ~polyMesh();


    // Member Functions

        // Database

            //- Over-ride the objectRegistry dbDir
            //  if this is single-region case
            virtual const fileName& dbDir() const;

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return the current instance directory for points
            //  Used in the consruction of gemometric mesh data dependent
            //  on points
            const word& pointsInstance() const;

            //- Return the current instance directory for cells
            //  Used in the consruction of tolpological mesh data dependent
            //  on cells
            const word& cellsInstance() const;


        // Access

            //- Return raw points
            const pointField& allPoints() const;

            //- Return raw faces
            const faceList& allFaces() const;

            //- Return raw cells
            const cellList& allCells() const;

            //- Return old mesh motion points
            const pointField& oldAllPoints() const;

            //- Return face owner
            const labelList& allOwner() const;

            //- Return face neighbour
            const labelList& allNeighbour() const;


            //- Return boundary mesh
            const polyBoundaryMesh& boundaryMesh() const
            {
                return boundary_;
            }

            //- Return point zone mesh
            const pointZoneMesh& pointZones() const
            {
                return pointZones_;
            }

            //- Return face zone mesh
            const faceZoneMesh& faceZones() const
            {
                return faceZones_;
            }

            //- Return cell zone mesh
            const cellZoneMesh& cellZones() const
            {
                return cellZones_;
            }

            //- Return parallel info
            const parallelInfo& parallelData() const;

            //- Return morph engine
            const polyMeshMorphEngine& morphEngine() const;


        // Constructor helper functions

            //- Add boundary patches
            void addPatches(const List<polyPatch*>&);

            //- Add mesh zones
            void addZones
            (
                const List<pointZone*>& pz,
                const List<faceZone*>& fz,
                const List<cellZone*>& cz
            );

            //- Add mesh modifiers to the morph engine
            void addTopologyModifiers(const List<polyMeshModifier*>& tm);

            //- Update the mesh based on the mesh files saved in
            //  time directories
            virtual readUpdateState readUpdate();


        // Mesh motion

            //- Is mesh moving
            bool moving() const
            {
                return moving_;
            }

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const pointField&);

            //- Reset motion
            void resetMotion() const;


        // Topological change

            //- Is mesh morphing topology
            bool morphing() const;

            //- Update mesh topology using the morph engine
            virtual void updateTopology();

            //- Return morph mapping data
            const mapPolyMesh& morphMap() const;

            //- Reset morphing
            void resetMorph() const;

            //- Return true if given face label is active in the mesh
            bool isActiveFace(const label faceIndex) const;


            // Special mesh functionality

                //- Purge polygonal face descriptions from extra points.
                //  keeping the points in point zones
                void purgePolygons();

                //- Correct face directions
                void correctFaceDirections();

                //- Zip up cells.  This function modifies the list of faces
                // such that all the cells are topologically closed.
                // Note that, while this operation is perfectly
                // sensible for correct meshes, it may lead to
                // arbitryary errors if the mesh is originally broken.
                // Use with care.  
                void zipUpCells();


        //  Storage management

            //- Print a list of all the currently allocated mesh data
            void printAllocated() const;

            //- Clear geometry
            void clearGeom();

            //- Clear addressing
            void clearAddressing();

            //- Clear all geometry and addressing unnecessary for CFD
            void clearOut();

            //- Clear primitive data (points, faces and cells)
            void clearPrimitives();

            //- Remove all files
            void removeFiles(const fileName& dir) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
