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
    blended

Description
    linear/upwind blended differencing scheme.

SourceFiles
    blendedMake.C

\*---------------------------------------------------------------------------*/

#ifndef blended_H
#define blended_H

#include "linear.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blended
:
    public linear<Type>,
    public upwind<Type>
{
    // Private data

        const scalar blendingFactor_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        blended(const blended&);

        //- Disallow default bitwise assignment
        void operator=(const blended&);


public:

    //- Runtime type information
    TypeName("blended");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        blended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const scalar blendingFactor
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            linear<Type>(mesh),
            upwind<Type>(mesh, faceFlux),
            blendingFactor_(blendingFactor)
        {}

        //- Construct from mesh and Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        blended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            linear<Type>(mesh),
            upwind<Type>
            (
                mesh,
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            ),
            blendingFactor_(readScalar(is))
        {}

        //- Construct from mesh, faceFlux and Istream
        blended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            linear<Type>(mesh),
            upwind<Type>(mesh, faceFlux),
            blendingFactor_(readScalar(is))
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return
                blendingFactor_*linear<Type>::weights(vf)
              + (1 - blendingFactor_)*upwind<Type>::weights(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
