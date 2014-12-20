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
    geometricTetPointFieldReconstructor

Description
    Tet point field reconstructor.

SourceFiles
    geometricTetPointFieldReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef geometricTetPointFieldReconstructor_H
#define geometricTetPointFieldReconstructor_H

#include "PtrList.H"
#include "Time.H"
#include "tetPolyMesh.H"
#include "tetPointMesh.H"
#include "IOobjectList.H"
#include "tetPolyPatchField.H"
#include "TetPointPatchFieldMapper.H"
#include "GeometricField.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class geometricTetPointFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class geometricTetPointFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        const tetPolyMesh& mesh_;

        //- List of processor meshes
        const PtrList<tetPolyMesh>& procMeshes_;

        //- List of processor point addressing lists
        const PtrList<labelIOList>& pointProcAddressing_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const PtrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        geometricTetPointFieldReconstructor
        (
            const geometricTetPointFieldReconstructor&
        );

        //- Disallow default bitwise assignment
        void operator=(const geometricTetPointFieldReconstructor&);

        //- Processor-to-global addressing for all points
        labelList procAddressing(const label procNo) const;

        //- Patch-to-patch addressing
        labelList procPatchAddressing
        (
            const labelList& procToGlobalAddr,
            const label procNo,
            const label patchNo
        ) const;


public:

        class tetPolyPatchFieldReconstructor
        :
            public TetPointPatchFieldMapper
        {
            label size_;

        public:

            // Constructors

                //- Construct given size
                tetPolyPatchFieldReconstructor(const label size)
                :
                    size_(size)
                {}


            // Member functions

                label size() const
                {
                    return size_;
                }

                bool direct() const
                {
                    return true;
                }

                const unallocLabelList& directAddressing() const
                {
                    return unallocLabelList::null();
                }
        };


    // Constructors

        //- Construct from components
        geometricTetPointFieldReconstructor
        (
            tetPolyMesh& mesh,
            const PtrList<tetPolyMesh>& procMeshes,
            const PtrList<labelIOList>& pointProcAddressing,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& cellProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Reconstruct field
        template<class Type>
        tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
        reconstructTetPointField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct and write all fields
        template<class Type>
        void reconstructTetPointFields
        (
            const IOobjectList& objects
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "geometricTetPointFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
