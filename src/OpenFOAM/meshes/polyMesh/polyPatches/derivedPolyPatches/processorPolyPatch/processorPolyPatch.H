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
    processorPolyPatch

Description
    Neighbour processor patch.

    Note: morph patch face ordering comes geometric or topological.
    Geometric: no cyclics allowed (assumes faces coincident)
    Topological: needs unmodified faces on both sides to correspond. Also
    needs at least one per connected patch area (so all patch faces can be
    visited from an unmodified face)

SourceFiles
    processorPolyPatch.C
    processorPolyPatchMorph.C

\*---------------------------------------------------------------------------*/

#ifndef processorPolyPatch_H
#define processorPolyPatch_H

#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class processorPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class processorPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        int myProcNo_;
        int neighbProcNo_;

        //- Processor-neighbbour patch face centres
        vectorField neighbFaceCentres_;

        //- Processor-neighbbour patch face areas
        vectorField neighbFaceAreas_;

        //- Processor-neighbbour patch neighbour cell centres
        vectorField neighbFaceCellCentres_;


    // Private static data

        //- Whether to use geometric or topological matching
        static bool geometricMatch_;

        //- Relative tolerance (for geometric matching only). Is factor of
        //  maximum edge length per face.
        static scalar matchTol_;


    // Private member functions

        // Morphing helper functions

            //- Do sending for geometric matching
            void sendGeometricOrder
            (
                const polyTopoChange&,
                const mapPolyMesh&
            ) const;

            //- Do ordering for geometric matching
            bool geometricOrder
            (
                const polyTopoChange&,
                const mapPolyMesh&,
                labelList&,
                labelList&
            ) const;

            //- Do sending for topological matching
            void sendTopologicalOrder
            (
                const polyTopoChange&,
                const mapPolyMesh&
            ) const;

            //- Do ordering for topological matching
            bool topologicalOrder
            (
                const polyTopoChange&,
                const mapPolyMesh&,
                labelList&,
                labelList&
            ) const;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initGeometry();

        //- Calculate the patch geometry
        void calcGeometry();

        //- Initialise the patches for moving points
        void initMovePoints(const pointField&);

        //- Correct patches after moving points
        void movePoints(const pointField&);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from components
        processorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo
        );

        //- Construct from Istream
        processorPolyPatch
        (
            Istream& is,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct from dictionary
        processorPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        processorPolyPatch(const processorPolyPatch&, const polyBoundaryMesh&);

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        processorPolyPatch
        (
            const processorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new processorPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new processorPolyPatch
                (
                    refCast<const processorPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~processorPolyPatch();


    // Member functions

        //- Return processor number
        int myProcNo() const
        {
            return myProcNo_;
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return neighbProcNo_;
        }

        //- Does the processor own the patch ?
        bool owner() const
        {
            return (myProcNo_ < neighbProcNo_);
        }

        //- Is the processor the patch neighbour ?
        bool neighbour() const
        {
            return !owner();
        }

        //- Force calculation of transformation tensors
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr
        ) const
        {
            coupledPolyPatch::calcTransformTensors(Cf, Cr, nf, nr);
        }

        //- Return processor-neighbbour patch face centres
        const vectorField& neighbFaceCentres() const
        {
            return neighbFaceCentres_;
        }

        //- Return processor-neighbbour patch face areas
        const vectorField& neighbFaceAreas() const
        {
            return neighbFaceAreas_;
        }

        //- Return processor-neighbbour patch neighbour cell centres
        const vectorField& neighbFaceCellCentres() const
        {
            return neighbFaceCellCentres_;
        }


        // Morph

            //- Initialize ordering (on new mesh)
            //  Gets -changes (in old mesh ordering)
            //  -map from old to new mesh
            virtual void sendOrder(const polyTopoChange&, const mapPolyMesh&)
             const;

            //- Morph:Return new ordering. Ordering is -faceMap: for every face
            //  index of the new face -rotation:for every new face the clockwise
            //  shift of the original face. Return false if nothing changes
            //  (faceMap is identity, rotation is 0), true otherwise.
            virtual bool order
            (
                const polyTopoChange&,
                const mapPolyMesh&,
                labelList& faceMap,
                labelList& rotation
            ) const;


        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
