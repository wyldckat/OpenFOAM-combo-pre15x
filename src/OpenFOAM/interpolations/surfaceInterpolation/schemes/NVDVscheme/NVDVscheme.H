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
    NVDVscheme

Description

    Class to create the weighting-factors based on the NVDV
    (Normalised Variable Diagram).
    The particular differencing scheme class is supplied as a template argument,
    the weight function of which is called by the weight function of this class
    for the internal faces as well as faces of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the face-flux, the cell and face
    gradients (from which the normalised variable distribution may be created)
    and the cell centre distance.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

SourceFiles
    NVDVscheme.C

\*---------------------------------------------------------------------------*/

#ifndef NVDVscheme_H
#define NVDVscheme_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class NVDVscheme Declaration
\*---------------------------------------------------------------------------*/

template<class NVDVweight>
class NVDVscheme
:
    public surfaceInterpolationScheme<vector>,
    public NVDVweight
{
    // Private data

        const surfaceScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        NVDVscheme(const NVDVscheme&);

        //- Disallow default bitwise assignment
        void operator=(const NVDVscheme&);


public:

    //- Define a typedef for the NVDVweight
    typedef NVDVweight Weight;


    //- Runtime type information
    TypeName("NVDVscheme");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        NVDVscheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const NVDVweight& weight
        )
        :
            surfaceInterpolationScheme<vector>(mesh),
            NVDVweight(weight),
            faceFlux_(faceFlux)
        {}

        //- Construct from mesh and Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        NVDVscheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<vector>(mesh),
            NVDVweight(is),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        NVDVscheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<vector>(mesh),
            NVDVweight(is),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<vector, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeNVDVsurfaceInterpolationScheme(SS, WEIGHT, NAME)                   \
                                                                               \
typedef NVDVscheme<WEIGHT> NVDVscheme##WEIGHT_;                                \
defineTemplateTypeNameAndDebugWithName(NVDVscheme##WEIGHT_, NAME, 0);          \
                                                                               \
surfaceInterpolationScheme<vector>::addMeshConstructorToTable                  \
<NVDVscheme<WEIGHT> > add##SS##vectorMeshConstructorToTable_;                  \
                                                                               \
surfaceInterpolationScheme<vector>::addMeshFluxConstructorToTable              \
<NVDVscheme<WEIGHT> > add##SS##vectorMeshFluxConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NVDVscheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
