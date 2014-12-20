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
    skewLinearEdgeInterpolation

Description
    Skewness-corrected interpolation scheme class derived from linear and
    returns linear weighting factors but also applies an explicit correction.

SourceFiles
    skewLinearEdgeInterpolationMake.C

\*---------------------------------------------------------------------------*/

#ifndef skewLinearEdgeInterpolation_H
#define skewLinearEdgeInterpolation_H

#include "linearEdgeInterpolation.H"
#include "gaussFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class skewLinearEdgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewLinearEdgeInterpolation
:
    public linearEdgeInterpolation<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewLinearEdgeInterpolation(const skewLinearEdgeInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const skewLinearEdgeInterpolation&);


public:

    //- Runtime type information
    TypeName("skewLinear");


    // Constructors

        //- Construct from mesh
        skewLinearEdgeInterpolation(const faMesh& mesh)
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh)
        {}


        //- Construct from mesh and Istream
        skewLinearEdgeInterpolation
        (
            const faMesh& mesh,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh)
        {}


        //- Construct from mesh, faceFlux and Istream
        skewLinearEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField&,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return this->mesh().skew();
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, faPatchField, edgeMesh> >
        correction
        (
            const GeometricField<Type, faPatchField, areaMesh>& vf
        ) const
        {
            const faMesh& mesh = this->mesh();

            tmp<GeometricField<Type, faPatchField, edgeMesh> > tsfCorr
            (
                new GeometricField<Type, faPatchField, edgeMesh>
                (
                    IOobject
                    (
                        vf.name(),
                        mesh.time().timeName(),
                        mesh()
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
                    mesh.skewCorrectionVectors()
                  & linearEdgeInterpolation
                    <
                        typename outerProduct
                        <
                            vector,
                            typename pTraits<Type>::cmptType
                        >::type
                    >
                    (mesh).interpolate
                    (
                        fa::gaussGrad<typename pTraits<Type>::cmptType>
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
