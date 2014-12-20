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
    fvPatchFieldMapperPatchRef

Description

\*---------------------------------------------------------------------------*/

#ifndef fvPatchFieldMapperPatchRef_H
#define fvPatchFieldMapperPatchRef_H

#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fvPatchFieldMapperPatchRef Declaration
\*---------------------------------------------------------------------------*/

class fvPatchFieldMapperPatchRef
:
    public fvPatchFieldMapper
{
    // Private data

        const fvPatch& sourcePatch_;
        const fvPatch& targetPatch_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvPatchFieldMapperPatchRef(const fvPatchFieldMapperPatchRef&);

        //- Disallow default bitwise assignment
        void operator=(const fvPatchFieldMapperPatchRef&);


public:

    // Constructors

        //- Construct from components
        fvPatchFieldMapperPatchRef(const fvPatch& source, const fvPatch& target)
        :
            sourcePatch_(source),
            targetPatch_(target)
        {}


    // Member functions

        const fvPatch& sourcePatch() const
        {
            return sourcePatch_;
        }

        const fvPatch& targetPatch() const
        {
            return targetPatch_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
