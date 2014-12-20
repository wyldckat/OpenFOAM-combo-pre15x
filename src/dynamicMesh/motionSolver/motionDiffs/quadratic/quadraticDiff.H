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
    quadraticDiff

Description
    Quadratic distance-based motion diffusion.

SourceFiles
    quadraticDiff.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticDiff_H
#define quadraticDiff_H

#include "linearDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class quadraticDiff Declaration
\*---------------------------------------------------------------------------*/

class quadraticDiff
:
    public linearDiff
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        quadraticDiff(const quadraticDiff&);

        //- Disallow default bitwise assignment
        void operator=(const quadraticDiff&);


public:

    //- Runtime type information
    TypeName("quadratic");


    // Constructors

        //- Construct from components
        quadraticDiff
        (
            const tetPolyMesh& tetMesh,
            const dictionary&
        );


    // Destructor

        virtual ~quadraticDiff();


    // Member Functions

        //- Return diffusivity field
        virtual tmp<elementScalarField> gamma() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
