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
    PhiScheme

Description

    Class to create the weighting-factors based on the face-flux.
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
    PhiScheme.C

\*---------------------------------------------------------------------------*/

#ifndef PhiScheme_H
#define PhiScheme_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PhiScheme Declaration
\*---------------------------------------------------------------------------*/

template<class PhiWeight>
class PhiScheme
:
    public surfaceInterpolationScheme<vector>,
    public PhiWeight
{
    // Private data

        const surfaceScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PhiScheme(const PhiScheme&);

        //- Disallow default bitwise assignment
        void operator=(const PhiScheme&);


public:

    //- Define a typedef for the PhiWeight
    typedef PhiWeight Weight;


    //- Runtime type information
    TypeName("PhiScheme");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        PhiScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const PhiWeight& weight
        )
        :
            surfaceInterpolationScheme<vector>(mesh),
            PhiWeight(weight),
            faceFlux_(faceFlux)
        {}

        //- Construct from mesh and Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        PhiScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<vector>(mesh),
            PhiWeight(is),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        PhiScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<vector>(mesh),
            PhiWeight(is),
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

#define makePhiSurfaceInterpolationScheme(SS, WEIGHT, NAME)                    \
                                                                               \
typedef PhiScheme<WEIGHT> Phischeme##WEIGHT_;                                  \
defineTemplateTypeNameAndDebugWithName(Phischeme##WEIGHT_, NAME, 0);           \
                                                                               \
surfaceInterpolationScheme<vector>::addMeshConstructorToTable                  \
<PhiScheme<WEIGHT> > add##SS##vectorMeshConstructorToTable_;                   \
                                                                               \
surfaceInterpolationScheme<vector>::addMeshFluxConstructorToTable              \
<PhiScheme<WEIGHT> > add##SS##vectorMeshFluxConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PhiScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
