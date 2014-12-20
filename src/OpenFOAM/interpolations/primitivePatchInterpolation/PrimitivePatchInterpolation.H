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
    PrimitivePatchInterpolation

Description
    Interpolation class within a primitive patch. Allows interpolation from
    points to faces and vice versa

SourceFiles
    PrimitivePatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef PrimitivePatchInterpolation_H
#define PrimitivePatchInterpolation_H

#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class PrimitivePatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Patch>
class PrimitivePatchInterpolation
{
    // Private data

        //- Reference to patch
        const Patch& patch_;


    // Demand-driven data

        //- Face-to-edge weights
        mutable scalarField* faceToEdgeWeightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PrimitivePatchInterpolation(const PrimitivePatchInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const PrimitivePatchInterpolation&);


        //- Face-to-edge weights
        const scalarField& faceToEdgeWeights() const;

        //- Make face-to-edge weights
        void makeFaceToEdgeWeights() const;

        //- Clear weights
        void clearWeights();


public:

    // Constructors

        //- Construct from PrimitivePatch
        PrimitivePatchInterpolation
        (
            const Patch& p
        )
        :
            patch_(p),
            faceToEdgeWeightsPtr_(NULL)
        {}


    // Destructor

        ~PrimitivePatchInterpolation();


    // Member Functions

        //- Interpolate from faces to points
        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Interpolate from points to faces
        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const Field<Type>& pf
        ) const;

        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const tmp<Field<Type> >& tpf
        ) const;

        //- Interpolate from faces to edges
        template<class Type>
        tmp<Field<Type> > faceToEdgeInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToEdgeInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PrimitivePatchInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
