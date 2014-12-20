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
    parallelInfo

Description
    Various mesh related information for a parallel run. Upon construction
    constructs all info by using parallel communication. 

    Requires: 
    - all processor patches to have correct ordering.
    - all processorPatches to have their transforms set ('makeTransforms')

    The shared point addressing is quite interesting. It gives on each processor
    the vertices that cannot be set using a normal swap on processor patches.
    These are the vertices that are shared between more than 2 processors.

    There is an issue with these shared vertices if they originate from
    cyclics (i.e. are now separated processor patches). They will all be
    mapped to the same global point (so even though the processor points are
    not on the same location) since topologically they are one and the same.

    So if you ask for sharedPoints() you get only one of the coordinates of
    the topologically shared points.

    All the hard work of these shared points is done by the globalPoints class.


SourceFiles
    parallelInfo.C
    parallelInfoMorph.C

\*---------------------------------------------------------------------------*/

#ifndef parallelInfo_H
#define parallelInfo_H

#include "polyMesh.H"
#include "Switch.H"
#include "processorTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class parallelInfo Declaration
\*---------------------------------------------------------------------------*/

class parallelInfo
:
    public processorTopology
{

    // Private class

        // To combineReduce a pointField. Just appends all lists.
        template <class T>
        class plusEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                label n = x.size();

                x.setSize(x.size() + y.size());

                forAll(y, i)
                {
                    x[n++] = y[i];
                }
            }
        };



    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        // Data related to the complete mesh

            //- Are there cyclic-parallel faces
            //  (interpolation/FEM problem in foam2.2).  
            Switch cyclicParallel_;

            //- Bounding box of complete mesh
            boundBox bb_;

            //- Total number of points in the complete mesh
            label nTotalPoints_;

            //- Total number of faces in the complete mesh
            label nTotalFaces_;

            //- Total number of cells in the complete mesh
            label nTotalCells_;



        // Processor patch addressing (be careful if not running in parallel!)

            //- List of processor patch labels 
            //  (size of list = number of processor patches)
            labelList processorPatches_;

            //- List of indices into processorPatches_ for each patch.
            //  Index = -1 for non-processor patches.
            //  (size of list = number of patches)
            labelList processorPatchIndices_;

            //- processorPatchIndices_ of the neighbours processor patches
            labelList processorPatchNeighbours_;


        // Globally shared point addressing

            //- Total number of global points
            label nGlobalPoints_;

            //- Indices of local points that are globally shared
            labelList sharedPointLabels_;

            //- Indices of globally shared points in the master list
            //  This list contains all the shared points in the mesh
            labelList sharedPointAddr_;

            //- Shared point global labels.
            //  Global point index for every local shared point.
            //  Only valid if constructed with this information or if
            //  pointProcAddressing read.
            mutable labelList* sharedPointGlobalLabelsPtr_;


    // Private Member Functions

        void initProcAddr();

        void calcSharedPoints() const;
        
        void calcGeomSharedPoints() const;
        



        //- Disallow default bitwise copy construct
        parallelInfo(const parallelInfo&);

        //- Disallow default bitwise assignment
        void operator=(const parallelInfo&);

public:

    //- Runtime type information
    ClassName("parallelInfo");


    // Static data members

        //- Geomtric tolerance (fraction of bounding box)
        static const Foam::scalar matchTol_;


    // Constructors

        //- Construct from mesh, derive rest (does parallel communication!)
        parallelInfo(const polyMesh& mesh);

        //- Old behaviour: read constructor given IOobject and a polyMesh
        //  reference. Only use this for testing!
        parallelInfo(const IOobject& io, const polyMesh& mesh);

        //- Construct from components. !To be removed when we trust above
        //  construct from mesh!
        parallelInfo
        (
            const polyMesh& mesh,
            const bool parallel,
            const bool cyclicParallel,
            const label nTotalPoints,
            const label nTotalFaces,
            const label nTotalCells,
            const label nGlobalPoints,
            const labelList& sharedPointLabels,
            const labelList& sharedPointAddr,
            const labelList& sharedPointGlobalLabels
        );


    // Destructor

        ~parallelInfo();

        //- Remove all demand driven data
        void clearOut();


    // Member Functions

        // Access

            //- Return the mesh reference
            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Does the mesh contain processor patches? (also valid when
            //  not running parallel)
            bool parallel() const
            {
                return processorPatches_.size() > 0;
            }

            const boundBox& bb() const
            {
                return bb_;
            }

            //- Does the mesh contain cyclic parallel faces (only valid when
            //  running parallel)
            bool cyclicParallel() const
            {
                return cyclicParallel_;
            }

            //- Return total number of points in decomposed mesh
            label nTotalPoints() const
            {
                return nTotalPoints_;
            }

            //- Return total number of faces in decomposed mesh
            label nTotalFaces() const
            {
                return nTotalFaces_;
            }

            //- Return total number of cells in decomposed mesh
            label nTotalCells() const
            {
                return nTotalCells_;
            }


        // Processor patch addressing (be careful when not running in parallel)

            //- Return list of processor patch labels 
            //  (size of list = number of processor patches)
            const labelList& processorPatches() const
            {
                return processorPatches_;
            }

            //- Return list of indices into processorPatches_ for each patch.
            //  Index = -1 for non-processor parches.
            //  (size of list = number of patches)
            const labelList& processorPatchIndices() const
            {
                return processorPatchIndices_;
            }

            //- Return processorPatchIndices of the neighbours
            //  processor patches. -1 if not running parallel.
            const labelList& processorPatchNeighbours() const
            {
                return processorPatchNeighbours_;
            }


        // Globally shared point addressing (only valid when running parallel)

            //- Return number of globally shared points
            label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }

            //- Return indices of local points that are globally shared
            const labelList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            //- Return addressing into the complete globally shared points
            //  list
            //  Note: It is assumed that a (never constructed) complete
            //  list of globally shared points exists.  The set of shared
            //  points on the current processor is a subset of all shared
            //  points. Shared point addressing gives the index in the
            //  list of all globally shared points for each of the locally
            //  shared points.
            const labelList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- Return shared point global labels. Tries to read
            //  'pointProcAddressing' and returns list or -1 if none
            //  available.
            const labelList& sharedPointGlobalLabels() const;

            //- Collect coordinates of shared points on all processors.
            //  (does parallel communication!)
            //  Note: not valid for cyclicParallel since shared cyclic points
            //  are merged into single global point. (use geometricSharedPoints
            //  instead)
            pointField sharedPoints() const;

            //- Like sharedPoints but keeps cyclic points separate.
            //  (does geometric merging; uses matchTol_*bb as merging tolerance)
            //  Use sharedPoints() instead.
            pointField geometricSharedPoints() const;


        // Edit

            //- Update for moving points.
            void movePoints(const pointField& newPoints);

            //- Change global mesh data given a topological change. Does a
            //  full parallel analysis to determine shared points and
            //  boundaries.
            void updateTopology(const mapPolyMesh&);


        // Write

            bool write() const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const parallelInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
