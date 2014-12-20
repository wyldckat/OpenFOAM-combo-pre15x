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
    motionSmoother

Description
    Given a displacement moves the mesh by scaling the displacement back
    until there are no more mesh errors. Holds displacement field
    (read upon construction since need boundary conditions) and scaling factor
    and optional patch number on which to scale back displacement.

    E.g.
        
        // Construct iterative mesh mover.
        motionSmoother meshMover(mesh, labelList(1, patchI));

        // Set wanted displacement:
        meshMover.displacement() = ..

        for (label iter = 0; iter < maxIter; iter++)
        {
            if (meshMover.scaleMesh(true))
            {
                Info<< "Successfully moved mesh" << endl;
                return true;
            }
        }

SourceFiles
    motionSmoother.C
    motionSmootherTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef motionSmoother_H
#define motionSmoother_H

#include "pointFields.H"
#include "labelHashSet.H"
#include "PackedList.H"
#include "indirectPrimitivePatch.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class twoDPointCorrector;

/*---------------------------------------------------------------------------*\
                           Class motionSmoother Declaration
\*---------------------------------------------------------------------------*/

class motionSmoother
{
    // Private class

        //! Can be removed once in PstreamCombineReduce
        template <class T>
        class minEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                x = min(x, y);
            }
        };


    // Private data

        //- Reference to polyMesh. Non-const since we move mesh.
        polyMesh& mesh_;

        //- Reference to pointMesh
        pointMesh& pMesh_;

        //- Reference to face subset of all adaptPatchIDs
        indirectPrimitivePatch& pp_;

        //- Indices of fixedValue patches that we're allowed to modify the
        // displacement on.
        const labelList adaptPatchIDs_;


        // Smoothing parameters

        //- Amount by which to reduce displacement on points in error.
        const scalar reduction_;

        //- Number of smoothing iterations on scale_ to distribute error.
        const label nSmoothScale_;


        // Mesh check params

            //- Max non-ortho
            const scalar maxNonOrtho_;

            //- Minimum allowed pyramid volume for accepting a mesh.
            const scalar minVol_;

            //- Max allowed concave angle (in degrees, 0 is straight,
            //  <0 would be convex)
            const scalar maxConcave_;

            //- Min face area
            const scalar minArea_;


        // Internal data

        //- Displacement field
        pointVectorField displacement_;

        //- Scale factor for displacement
        pointScalarField scale_;

        //- Starting mesh position
        pointField oldPoints_;

        //- Is mesh point on boundary or not
        PackedList<1> isInternalPoint_;

        //- Is motion limited to twoD
        Switch twoDMotion_;

        //- Only if twoDMotion_ : 2-D motion corrector pointer
        twoDPointCorrector* correct2DPtr_;
        

    // Private Member Functions

        //- (unweighted) average of all values on points connected via edges
        //  to pointI
        template <class Type>
        static Type avg
        (
            const GeometricField<Type, pointPatchField, pointMesh>& fld,
            const edgeList& edges,
            const pointField& points,
            const labelList& edgeLabels,
            const label pointI
        );

        //- Distance weighted average. Is average with inverse distance
        //  weighting and varying edge-diffusivity.
        template <class Type>
        static Type avg
        (
            const GeometricField<Type, pointPatchField, pointMesh>& fld,
            const scalarField& edgeGamma,
            const edgeList& edges,
            const pointField& points,
            const labelList& edgeLabels,
            const label pointI
        );

        static void checkFld(const pointScalarField& fld);

        //- Get points used by given faces
        labelHashSet getPoints(const labelHashSet&) const;

        //- explicit smoothing and min on all internal points
        void minSmooth(const pointScalarField&, pointScalarField&) const;

        //- same but only on selected points (usually patch points)
        void minSmooth
        (
            const labelList& meshPoints,
            const pointScalarField& fld,
            pointScalarField& newFld
        ) const;

        //- Scale certain (internal) points of a field
        void scaleField
        (
            const labelHashSet& pointLabels,
            const scalar scale,
            pointScalarField& fld
        ) const;

        //- As above but points have to be in meshPoints as well
        //  (usually to scale patch points)
        void scaleField
        (
            const labelList& meshPoints,
            const labelHashSet& pointLabels,
            const scalar scale,
            pointScalarField& fld
        ) const;

        //- Helper function. Is point internal?
        bool isInternalPoint(const label pointI) const;

        //- Disallow default bitwise copy construct
        motionSmoother(const motionSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const motionSmoother&);


public:

    ClassName("motionSmoother");

    // Constructors

        //- Construct from mesh, optional fixedValue patch
        //  where displacement will be adapted and parameters: number of
        //  error distribution cycles, minimum pyramid vol, max. concaveness
        //  check. Will read displacement field.
        motionSmoother
        (
            polyMesh&,
            pointMesh&,
            indirectPrimitivePatch& pp,         // 'outside' points
            const labelList& adaptPatchIDs,     // patches forming 'outside'
            const scalar reduction = 0.5,       // scale back displacement
            const label nSmoothScale = 1,       // number of smoothing iters
            const scalar maxNonOrtho = 90,      // max non-ortho
            const scalar minVol = SMALL,        // min allowed cell vol
            const scalar maxConcave = 180,      // max allowed face concavity
            const scalar minArea = VSMALL       // min face area
        );


    // Destructor

        ~motionSmoother();


    // Member Functions

        // Access

            //- Reference to mesh
            const polyMesh& mesh() const;

            //- Reference to pointMesh
            const pointMesh& pMesh() const;

            //- Reference to patch
            const indirectPrimitivePatch& patch() const;

            scalar maxNonOrtho() const;

            scalar minVol() const;

            scalar maxConcave() const;

            scalar minArea() const;


            //- Reference to displacement field
            pointVectorField& displacement();

            //- Reference to displacement field
            const pointVectorField& displacement() const;

            //- Reference to scale field
            const pointScalarField& scale() const;

            //- Starting mesh position
            const pointField& oldPoints() const;

            bool twoDMotion() const
            {
                return twoDMotion_;
            }

            twoDPointCorrector& twoDCorrector()
            {
                if (!correct2DPtr_)
                {
                    FatalErrorIn("motionSmoother::correct2D()")
                        << "No 2D motion. Check constant/motionProperties"
                        <<  " and empty 3D boundaries" << abort(FatalError);
                }
                return *correct2DPtr_;
            }



        // Edit

            //- Take over existing mesh position.
            void correct();


            //- Move mesh. Does 2D correction (modifies passed pointField) and
            //  polyMesh::movePoints. Returns swept volumes.
            tmp<scalarField> movePoints(pointField&);

            //- Move mesh with given scale. Return true if mesh ok or has
            //  less than nAllow errors, false
            //  otherwise and locally update scale. Smoothmesh=false means only
            //  patch points get moved.
            //  Parallel ok (as long as displacement field is consistent
            //  across patches)
            bool scaleMesh
            (
                const bool smoothMesh = true,
                const label nAllow = 0
            );

            //- Update topology
            void updateTopology();

            //- Check mesh with current settings. Collects incorrect faces
            //  in set. Returns true if one or more faces in error.
            //  Parallel ok.
            bool checkMesh(labelHashSet&) const;


            // Helper functions to manipulate displacement vector.

                //- Point-jacobi smoothing of internal points
                template <class Type>
                void smooth
                (
                    GeometricField<Type, pointPatchField, pointMesh>&
                ) const;

                //- Fully explicit smoothing of internal points with varying
                //  diffusivity.
                template <class Type>
                void smooth
                (
                    const GeometricField<Type, pointPatchField, pointMesh>& fld,
                    const scalarField& edgeGamma,
                    GeometricField<Type, pointPatchField, pointMesh>& newFld
                ) const;

                //- Sychronizes patch points on pointField
                template<class Type, class CombineOp>
                static void syncField
                (
                    GeometricField<Type, pointPatchField, pointMesh>&,
                    const Type& zero,
                    const CombineOp& cop
                );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "motionSmootherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
