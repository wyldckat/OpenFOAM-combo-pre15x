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
    cubicCorrection

Description
    Cubic interpolation scheme class derived from linear and returns
    linear weighting factors but also applies an explicit correction.

SourceFiles
    cubicCorrectionMake.C

\*---------------------------------------------------------------------------*/

#ifndef cubicCorrection_H
#define cubicCorrection_H

#include "linear.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cubicCorrection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cubicCorrection
:
    public linear<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        cubicCorrection(const cubicCorrection&);

        //- Disallow default bitwise assignment
        void operator=(const cubicCorrection&);


public:

    //- Runtime type information
    TypeName("cubicCorrection");


    // Constructors

        //- Construct from mesh
        cubicCorrection(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh),
            linear<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        cubicCorrection
        (
            const fvMesh& mesh,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            linear<Type>(mesh)
        {}

        //- Construct from mesh, faceFlux and Istream
        cubicCorrection
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            linear<Type>(mesh)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            // calculate the appropriate interpolation factors
            const surfaceScalarField& lambda = mesh.weights();
            surfaceScalarField kSc = lambda*(1.0 - lambda*(3.0 - 2.0*lambda));
            surfaceScalarField kVecP = lambda*(1.0 - lambda*(2.0 - lambda));
            surfaceScalarField kVecN = sqr(lambda)*(lambda - 1.0);

            tmp<GeometricField<Type, fvPatchField, surfaceMesh> > tsfCorr
            (
                new GeometricField<Type, fvPatchField, surfaceMesh>
                (
                    IOobject
                    (
                        vf.name(),
                        mesh.time().timeName(),
                        mesh
                    ),
                    surfaceInterpolationScheme<Type>::interpolate(vf, kSc, -kSc)
                )
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                tsfCorr().replace
                (
                    cmpt,
                    tsfCorr().component(cmpt)
                  + (
                        linear
                        <
                            typename outerProduct
                            <
                                vector,
                                typename pTraits<Type>::cmptType
                            >::type
                        >(mesh).interpolate
                        (
                            fv::gaussGrad
                            <typename pTraits<Type>::cmptType>(mesh)
                           .grad(vf.component(cmpt)),
                            kVecP,
                            kVecN
                        ) & mesh.Sf()
                    )/mesh.magSf()/mesh.surfaceInterpolation::deltaCoeffs()
                );
            }

            return tsfCorr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
