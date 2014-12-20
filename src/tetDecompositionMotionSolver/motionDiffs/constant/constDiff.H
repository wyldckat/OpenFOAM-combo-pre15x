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
    constDiff

Description
    Constant motion diffusion.

SourceFiles
    constDiff.C

\*---------------------------------------------------------------------------*/

#ifndef constDiff_H
#define constDiff_H

#include "motionDiff.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class constDiff Declaration
\*---------------------------------------------------------------------------*/

class constDiff
:
    public motionDiff
{
private:

    // Private data
    
        scalar value_;
    
        elementScalarField gamma_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constDiff(const constDiff&);

        //- Disallow default bitwise assignment
        void operator=(const constDiff&);


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constDiff
        (
            const tetDecompositionMotionSolver& mSolver
        );


    // Destructor

        virtual ~constDiff();


    // Member Functions

        //- Return diffusivity field
        const elementScalarField& gamma() const
        {
            return gamma_;
        }

        //- Correct the motion diffusivity
        void correct()
	{}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
