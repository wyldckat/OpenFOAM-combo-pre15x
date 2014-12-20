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
    skewLinear

Description
    Skewness-corrected interpolation scheme class derived from linear and
    returns linear weighting factors but also applies an explicit correction.

SourceFiles
    skewLinear.C

\*---------------------------------------------------------------------------*/

#ifndef skewLinear_H
#define skewLinear_H

#include "skewCorrectionVectors.H"
#include "linear.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class skewLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewLinear
:
    public linear<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewLinear(const skewLinear&);

        //- Disallow default bitwise assignment
        void operator=(const skewLinear&);


public:

    //- Runtime type information
    TypeName("skewLinear");


    // Constructors

        //- Construct from mesh
        skewLinear(const fvMesh& mesh)
        :
            linear<Type>(mesh)
        {}


        //- Construct from mesh and Istream
        skewLinear
        (
            const fvMesh& mesh,
            Istream&
        )
        :
            linear<Type>(mesh)
        {}


        //- Construct from mesh, faceFlux and Istream
        skewLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            linear<Type>(mesh)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return skewCorrectionVectors::New(this->mesh()).skew();
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh);

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
                    mesh,
                    dimensioned<Type>
                    (
                        vf.name(),
                        vf.dimensions(),
                        pTraits<Type>::zero
                    )
                )
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                tsfCorr().replace
                (
                    cmpt,
                    scv() & linear
                    <
                        typename outerProduct
                        <
                            vector,
                            typename pTraits<Type>::cmptType
                        >::type
                    > (mesh).interpolate
                    (
                        fv::gaussGrad<typename pTraits<Type>::cmptType>
                        (mesh).grad(vf.component(cmpt))
                    )
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
