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
    inverseLinear

Description
    Inversed weight central-differencing interpolation scheme class.
    Useful for inverse weighted and harmonic interpolations.

SourceFiles
    inverseLinear.C

\*---------------------------------------------------------------------------*/

#ifndef inverseLinear_H
#define inverseLinear_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inverseLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class inverseLinear
:
    virtual public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const inverseLinear&);


public:

    //- Runtime type information
    TypeName("inverseLinear");


    // Constructors

        //- Construct from mesh
        inverseLinear(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        inverseLinear(const fvMesh& mesh, Istream&)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        inverseLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            const fvMesh& mesh = this->mesh();

            tmp<surfaceScalarField> tcdWeights
            (
                mesh.surfaceInterpolation::weights()
            );
            surfaceScalarField& cdWeights = tcdWeights();

            cdWeights.internalField() = 1.0 - cdWeights.internalField();

            forAll (mesh.boundary(), patchI)
            {
                if (mesh.boundary()[patchI].coupled())
                {
                  cdWeights.boundaryField()[patchI] =
                      1.0 - cdWeights.boundaryField()[patchI];
                }
            }

            return tcdWeights;
        }
};


template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
inverseLinearInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
    return surfaceInterpolationScheme<Type>::interpolate
    (
        vf,
        1.0 - vf.mesh().surfaceInterpolation::weights()
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
inverseLinearInterpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, surfaceMesh> > tinterp = 
        inverseLinearInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
