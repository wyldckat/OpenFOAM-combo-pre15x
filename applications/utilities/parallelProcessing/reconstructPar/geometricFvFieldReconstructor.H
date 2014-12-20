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
    geometricFvFieldReconstructor

Description
    FV volume and surface field reconstructor.

SourceFiles
    geometricFvFieldReconstructor.C
    geometricFvFieldReconstructorReconstructFields.C

\*---------------------------------------------------------------------------*/

#ifndef geometricFvFieldReconstructor_H
#define geometricFvFieldReconstructor_H

#include "ptrList.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "fvPatchFieldMapper.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                 Class geometricFvFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class geometricFvFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        fvMesh& mesh_;

        //- List of processor meshes
        const ptrList<fvMesh>& procMeshes_;

        //- List of processor point addressing lists
        const ptrList<labelIOList>& pointProcAddressing_;

        //- List of processor face addressing lists
        const ptrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const ptrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const ptrList<labelIOList>& boundaryProcAddressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        geometricFvFieldReconstructor(const geometricFvFieldReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const geometricFvFieldReconstructor&);


public:

        class fvPatchFieldReconstructor
        :
            public fvPatchFieldMapper
        {
            labelList directAddressing_;

        public:

            // Constructors

                //- Construct given addressing
                fvPatchFieldReconstructor(const label size)
                :
                    directAddressing_(size, 0)
                {}


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


    // Constructors

        //- Construct from components
        geometricFvFieldReconstructor
        (
            fvMesh& mesh,
            const ptrList<fvMesh>& procMeshes,
            const ptrList<labelIOList>& pointProcAddressing,
            const ptrList<labelIOList>& faceProcAddressing,
            const ptrList<labelIOList>& cellProcAddressing,
            const ptrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Reconstruct volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        reconstructFvVolumeField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct surface field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
        reconstructFvSurfaceField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct and write all volume fields
        template<class Type>
        void reconstructFvVolumeFields
        (
            const IOobjectList& objects
        );

        //- Reconstruct and write all volume fields
        template<class Type>
        void reconstructFvSurfaceFields
        (
            const IOobjectList& objects
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "geometricFvFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
