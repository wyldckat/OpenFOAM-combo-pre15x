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
    patchEnhancedDiff

Description
    Patch-enhanced motion diffusion.

SourceFiles
    patchEnhancedDiff.C

\*---------------------------------------------------------------------------*/

#ifndef patchEnhancedDiff_H
#define patchEnhancedDiff_H

#include "motionDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class patchEnhancedDiff Declaration
\*---------------------------------------------------------------------------*/

class patchEnhancedDiff
:
    public motionDiff
{
private:

    // Private data

        //- Patches selected to base the distance on
        wordList patchNames_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchEnhancedDiff(const patchEnhancedDiff&);

        //- Disallow default bitwise assignment
        void operator=(const patchEnhancedDiff&);


public:

    //- Runtime type information
    TypeName("patchEnhanced");


    // Constructors

        //- Construct from components
        patchEnhancedDiff
        (
            const tetPolyMesh& tetMesh,
            const dictionary&
        );


    // Destructor

        virtual ~patchEnhancedDiff();


    // Member Functions

        //- Enhance given difusivity
        void enhance(elementScalarField&) const;

        //- Return diffusivity field
        virtual tmp<elementScalarField> gamma() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
