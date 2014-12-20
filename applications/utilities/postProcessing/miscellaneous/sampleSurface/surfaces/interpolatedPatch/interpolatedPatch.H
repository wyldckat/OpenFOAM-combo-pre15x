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
    interpolatedPatch

Description

SourceFiles
    interpolatedPatch.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedPatch_H
#define interpolatedPatch_H

#include "surface.H"
#include "constantPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshSearch;
class dictionary;
class plane;
class volPointInterpolation;
template<class Type> class fieldsCache;

/*---------------------------------------------------------------------------*\
                           Class interpolatedPatch Declaration
\*---------------------------------------------------------------------------*/

class interpolatedPatch
:
    public constantPatch
{
    // Private Member Functions

        //- Interpolate field onto faces
        template <class T>
        Field<T> doInterpolate
        (
            const word& fieldName,
            const fieldsCache<T>& cache,
            const volPointInterpolation& pInterp,
            const dictionary& interpolationSchemes
        ) const;

public:

    //- Runtime type information
    TypeName("interpolatedPatch");

    // Constructors

        //- Construct from components
        interpolatedPatch
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& name,
            const word& patchName,
            const bool triangulate = true
        );

        //- Construct from dictionary
        interpolatedPatch
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict          
        );


    // Destructor

        virtual ~interpolatedPatch();


    // Member Functions

        //- correct() implemented by constantPatch::correct()

        //- interpolate field to surface
        virtual scalarField interpolate
        (
            const word&,
            const fieldsCache<scalar>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual vectorField interpolate
        (
            const word&,
            const fieldsCache<vector>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tensorField interpolate
        (
            const word&,
            const fieldsCache<tensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
