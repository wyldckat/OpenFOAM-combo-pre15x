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
    Gamma

Description
    Class to create the weighting-factors based on the C-NVD
    (Courant number limited Normalised Variable Diagram).
    The particular differencing scheme class is supplied as a template argument,
    the weight function of which is called by the weight function of this class
    for the internal faces as well as faces of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the face-flux, the cell variable
    values, the cell gradients the face Courant number and the cell centre
    distance.  This version of the CNVD does not use the face gradient in the
    creation of the normalised variables but used the cell-centre variable
    values directly.  This allows a limit to be applied to the estimated
    upwind-cell values as required by the CICSAM scheme of Onno Ubbink.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

SourceFiles
    CNVDscheme.C

\*---------------------------------------------------------------------------*/

#ifndef CNVDscheme_H
#define CNVDscheme_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CNVDscheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class CNVDweight>
class CNVDscheme
:
    public surfaceInterpolationScheme<Type>,
    public CNVDweight
{
    // Private data

        const surfaceScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CNVDscheme(const CNVDscheme&);

        //- Disallow default bitwise assignment
        void operator=(const CNVDscheme&);


public:

    //- Define a typedef for the CNVDweight
    typedef CNVDweight Weight;


    //- Runtime type information
    TypeName("CNVDscheme");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        CNVDscheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const CNVDweight& weight
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            CNVDweight(weight),
            faceFlux_(faceFlux)
        {}

        //- Construct from mesh and Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        CNVDscheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            CNVDweight(is),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        CNVDscheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            CNVDweight(is),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeCNVDsurfaceInterpolationTypeScheme(SS, WEIGHT, NAME, TYPE)         \
                                                                               \
typedef CNVDscheme<TYPE, WEIGHT> CNVDscheme##TYPE##WEIGHT_;                    \
defineTemplateTypeNameAndDebugWithName(CNVDscheme##TYPE##WEIGHT_, NAME, 0);    \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<CNVDscheme<TYPE, WEIGHT> >                                                    \
    add##SS##TYPE##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<CNVDscheme<TYPE, WEIGHT> >                                                    \
    add##SS##TYPE##MeshFluxConstructorToTable_;

#define makeCNVDsurfaceInterpolationScheme(SS, WEIGHT, NAME)                   \
                                                                               \
makeCNVDsurfaceInterpolationTypeScheme(SS, WEIGHT, NAME, scalar)               \
makeCNVDsurfaceInterpolationTypeScheme(SS, WEIGHT, NAME, vector)               \
makeCNVDsurfaceInterpolationTypeScheme(SS, WEIGHT, NAME, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CNVDscheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
