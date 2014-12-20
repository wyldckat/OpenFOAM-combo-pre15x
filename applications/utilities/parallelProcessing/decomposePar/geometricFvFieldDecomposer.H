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
    geometricFvFieldDecomposer

Description
    Finite Volume volume and surface field decomposer.

SourceFiles
    geometricFvFieldDecomposer.C
    geometricFvFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef geometricFvFieldDecomposer_H
#define geometricFvFieldDecomposer_H

#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class geometricFvFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class geometricFvFieldDecomposer
{
    // Private data

        //- Reference to processor mesh
        const fvMesh& processorMesh_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        geometricFvFieldDecomposer(const geometricFvFieldDecomposer&);

        //- Disallow default bitwise assignment
        void operator=(const geometricFvFieldDecomposer&);


public:

        //- Patch field decomposer class
        class fvPatchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            labelList directAddressing_;

        public:

            // Constructors

                //- Construct given addressing
                fvPatchFieldDecomposer
                (
                    const unallocLabelList& addressingSlice,
                    const label addressingOffset
                )
                :
                    directAddressing_(addressingSlice)
                {
                    forAll (directAddressing_, i)
                    {
                        // Subtract one to align addressing.  
                        directAddressing_[i] -= addressingOffset + 1;
                    }
                }


            // Member functions

                label size() const
                {
                    return directAddressing_.size();
                }

                virtual bool resizeOnly() const
                {
                    return false;
                }

                bool direct() const
                {
                    return true;
                }

                const unallocLabelList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class
        class processorFvPatchVolFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            labelListList addressing_;
            scalarListList weights_;

        public:

            //- Construct given addressing
            processorFvPatchVolFieldDecomposer
            (
                const fvMesh& mesh,
                const unallocLabelList& addressingSlice
            )
            :
                addressing_(addressingSlice.size()),
                weights_(addressingSlice.size())
            {
                const scalarField& weights = mesh.weights().internalField();
                const labelList& own = mesh.faceOwner();
                const labelList& neighb = mesh.faceNeighbour();

                forAll (addressing_, i)
                {
                    // Subtract one to align addressing.  
                    label ai = mag(addressingSlice[i]) - 1;

                    if (ai < neighb.size())
                    {
                        // This is a regular face. it has been an internal face
                        // of the original mesh and now it has become a face
                        // on the parallel boundary
                        addressing_[i].setSize(2);
                        weights_[i].setSize(2);

                        addressing_[i][0] = own[ai];
                        addressing_[i][1] = neighb[ai];

                        weights_[i][0] = weights[ai];
                        weights_[i][1] = 1.0 - weights[ai];
                    }
                    else
                    {
                        // This is a face that used to be on a cyclic boundary
                        // but has now become a parallel patch face. I cannot
                        // do the interpolation properly (I would need to look
                        // up the different (face) list of data), so I will
                        // just grab the value from the owner cell
                        // 
                        addressing_[i].setSize(1);
                        weights_[i].setSize(1);

                        addressing_[i][0] = own[ai];

                        weights_[i][0] = 1.0;
                    }
                }
            }


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual bool resizeOnly() const
                {
                    return false;
                }

                bool direct() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


        //- Processor patch field decomposer class
        class processorFvPatchSurfaceFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            labelListList addressing_;
            scalarListList weights_;

        public:

            //- Construct given addressing
            processorFvPatchSurfaceFieldDecomposer
            (
                const unallocLabelList& addressingSlice
            )
            :
                addressing_(addressingSlice.size()),
                weights_(addressingSlice.size())
            {
                forAll (addressing_, i)
                {
                    addressing_[i].setSize(1);
                    weights_[i].setSize(1);

                    addressing_[i][0] = mag(addressingSlice[i]) - 1;
                    weights_[i][0] = sign(addressingSlice[i]);
                }
            }


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual bool resizeOnly() const
                {
                    return false;
                }

                bool direct() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


    // Constructors

        //- Construct from components
        geometricFvFieldDecomposer
        (
            const fvMesh& processorMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );


    // Member Functions

        // Read the fields and hold on the pointer list
        template<class GeoField>
        static void readFields
        (
            const fvMesh& mesh,
            const IOobjectList& objects,
            PtrList<GeoField>& fields
        );

        //- Decompose volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        decomposeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
        decomposeField
        (
            const GeometricField<Type, fvPatchField, surfaceMesh>& field
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "geometricFvFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
