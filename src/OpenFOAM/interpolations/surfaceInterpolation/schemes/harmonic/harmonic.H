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
    harmonic

Description
    Harmonic-mean differencing scheme class.  This scheme interpolates 1/field
    using a scheme specified at run-time and return the reciprocal of the
    interpolate.

SourceFiles
    harmonic.C

\*---------------------------------------------------------------------------*/

#ifndef harmonic_H
#define harmonic_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "inverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class harmonic Declaration
\*---------------------------------------------------------------------------*/

class harmonic
:
    virtual public surfaceInterpolationScheme<scalar>
{
    // Private data

        //- The basic linear scheme 
        tmp<surfaceInterpolationScheme<scalar> > tlinearScheme_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const harmonic&);


public:

    //- Runtime type information
    TypeName("harmonic");


    // Constructors

        //- Construct from Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        harmonic
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh),
            tlinearScheme_
            (
             // Hard-coding the interpolation scheme to be inverse linear
//                 surfaceInterpolationScheme<scalar>::New(mesh, is)
                new inverseLinear<scalar>(mesh, is)
            )
        {}

        //- Construct from faceFlux and Istream
        harmonic
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh),
            tlinearScheme_
            (
             // Hard-coding the interpolation scheme to be inverse linear
//                 surfaceInterpolationScheme<scalar>::New(mesh, faceFlux, is)
                new inverseLinear<scalar>(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<scalar, fvPatchField, volMesh>&
        ) const
        {
            notImplemented
            (
                "harmonic::weights"
                "(const GeometricField<scalar, fvPatchField, volMesh>&)"
            );

            return tmp<surfaceScalarField>(NULL);
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<scalar, fvPatchField, surfaceMesh> >
        interpolate
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vf
        ) const
        {
            return 1.0/(tlinearScheme_().interpolate(1.0/vf));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
