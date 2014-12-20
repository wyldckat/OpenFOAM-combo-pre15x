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
    meshWave

Description
    Wave propagation of information through grid. Every iteration
    information goes through one layer of cells. Templated on information
    that is transferred.

    Handles parallel and cyclics and non-parallel cyclics.

    Note: whether to propagate depends on the return value of Type::update
    which returns true (i.e. propagate) if the value changes by more than a
    certain tolerance.
    This tolerance can be very strict for normal face-cell and parallel
    cyclics (we use a value of 0.01 just to limit propagation of small changes)
    but for non-parallel cyclics this tolerance can be critical and if chosen
    too small can lead to non-convergence.

SourceFiles
    meshWave.C

\*---------------------------------------------------------------------------*/

#ifndef meshWave_H
#define meshWave_H

#include "label.H"
#include "boolList.H"
#include "labelList.H"
#include "scalarField.H"
#include "volFields.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class meshWaveName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(meshWave);


/*---------------------------------------------------------------------------*\
                           Class meshWave Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class meshWave
:
    public meshWaveName
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Wall information for all cells
        List<Type> allCellInfo_;

        //- Wall information for all faces
        List<Type> allFaceInfo_;

        //- Has face changed
        boolList changedFace_;

        //- List of changed faces
        labelList changedFaces_;

        //- Number of changed faces
        label nChangedFaces_;

        // Cells that have changed
        boolList changedCell_;
        labelList changedCells_;
        label nChangedCells_;

        //- Contains cyclics
        bool hasCyclicPatches_;

        //- Contains processor patches
        bool hasProcPatches_;

        //- Number of evaluations
        label nEvals_;

        //- Number of unvisited cells/faces
        label nUnvisitedCells_;
        label nUnvisitedFaces_;

        //- Iteration counter
        label iter_;

        //- For debugging: string containing processor number
        const string procLabel_;

    // Static Functions

        //- Write faces info
        static Ostream& writeFaces
        (
            const label nFaces,
            const labelList& faceLabels,
            const List<Type>& faceInfo,
            Ostream& os
        );

        //- Read faces info
        static Istream& readFaces
        (
            label& nFaces,
            labelList& faceLabels,
            List<Type>& faceInfo,
            Istream& is
        );


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshWave(const meshWave&);

        //- Disallow default bitwise assignment
        void operator=(const meshWave&);


        //- Updates cellInfo with information from neighbour. Updates all
        //  statistics.
        bool updateCell
        (
            const label cellI,
            const label neighbourFaceI,
            const Type& neighbourInfo,
            const scalar tol,
            Type& cellInfo
        );

        //- Updates faceInfo with information from neighbour. Updates all
        //  statistics.
        bool updateFace
        (
            const label faceI,
            const label neighbourCellI,
            const Type& neighbourInfo,
            const scalar tol,
            Type& faceInfo
        );

        //- Updates faceInfo with information from same face. Updates all
        //  statistics.
        bool updateFace
        (
            const label faceI,
            const Type& neighbourInfo,
            const scalar tol,
            Type& faceInfo
        );

        //- Copy face info into member data
        void setFaceInfo
        (
            const labelList& changedFaces,
            const List<Type>& changedFacesInfo
        );


        // Parallel, cyclic

            //- Debugging: check info on both sides of cyclic
            void checkCyclic(const polyPatch& pPatch) const;

            //- Has patches of certain type?
            bool hasPatchType(const word& nameOfType);

            //- Merge received patch data into global data
            void mergeFaceInfo
            (
                const polyPatch& patch,
                const label nFaces,
                const labelList&,
                const List<Type>&,
                const bool isParallel
            );

            //- Extract info for single patch only
            label getChangedPatchFaces
            (
                const polyPatch& patch,
                const label startFaceI,
                const label nFaces,
                labelList& changedPatchFaces,
                List<Type>& changedPatchFacesInfo
            ) const;

            //- Handle leaving domain. Implementation referred to Type
            void leaveDomain
            (
                const polyPatch& patch,
                const label nFaces,
                const labelList& faceLabels,
                List<Type>& faceInfo
            ) const;

            //- Handle leaving domain. Implementation referred to Type
            void enterDomain
            (
                const polyPatch& patch,
                const label nFaces,
                const labelList& faceLabels,
                List<Type>& faceInfo
            ) const;

            //- Send info to neighbour
            void sendPatchInfo
            (
                const label neighbour,
                const label nFaces,
                const labelList&,
                const List<Type>&
            ) const;

            //- Receive info from neighbour. Returns number of faces received.
            label receivePatchInfo
            (
                const label neighbour,
                labelList&,
                List<Type>&
            ) const;

            //- Offset face labels by constant value
            static void offset
            (
                const polyPatch& patch,
                const label off,
                const label nFaces,
                labelList& faces
            );

            //- Apply transformation to Type
            static void transform
            (
                const tensorField& rotTensor,
                const label nFaces,
                List<Type>& faceInfo
            );

            //- Merge data from across processor boundaries
            void handleProcPatches();

            //- Merge data from across cyclics
            void handleCyclicPatches();


      // Private static data

            static const scalar geomTol_;
            static const scalar propagationTol_;

public:

    // Constructors

        //- Construct from mesh and list of changed faces with the Type
        //  for these faces. Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        meshWave
        (
            const polyMesh& mesh,
            const labelList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            const label maxIter
        );

        //- Construct from mesh, list of changed faces with the Type
        //  for these faces and initial field.
        //  Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        meshWave
        (
            const polyMesh& mesh,
            const labelList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            const List<Type>& allCellInfo,
            const label maxIter
        );


    // Destructor

        ~meshWave();


    // Member Functions

        //- Get allFaceInfo
        const List<Type>& allFaceInfo() const
        {
            return allFaceInfo_;
        }

        //- Get allCellInfo
        const List<Type>& allCellInfo() const
        {
            return allCellInfo_;
        }

        //- Get number of unvisited cells, i.e. cells that were not (yet)
        //  reached from walking across mesh. This can happen from
        //  - not enough iterations done
        //  - a disconnected mesh
        //  - a mesh without walls in it
        label getUnsetCells() const;

        label getUnsetFaces() const;

        //- Propagate from face to cell. Returns total number of cells
        //  (over all processors) changed.
        label faceToCell();

        //- Propagate from cell to face. Returns total number of faces
        //  (over all processors) changed. (Faces on processorpatches are
        //  counted double)
        label cellToFace();

        //- Iterate until no changes or maxIter reached.
        label iterate(const label maxIter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
