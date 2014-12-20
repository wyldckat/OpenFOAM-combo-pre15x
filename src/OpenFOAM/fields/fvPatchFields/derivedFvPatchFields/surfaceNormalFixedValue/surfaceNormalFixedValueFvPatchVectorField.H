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
    surfaceNormalFixedValueFvPatchVectorField

Description
    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    surfaceNormalFixedValueFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceNormalFixedValueFvPatchVectorField_H
#define surfaceNormalFixedValueFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class surfaceNormalFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class surfaceNormalFixedValueFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    
        scalarField refValue_;


public:

    //- Runtime type information
    TypeName("surfaceNormalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        surfaceNormalFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&
        );

        //- Construct from patch, internal field and dictionary
        surfaceNormalFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceNormalFixedValueFvPatchVectorField
        //  onto a new patch
        surfaceNormalFixedValueFvPatchVectorField
        (
            const surfaceNormalFixedValueFvPatchVectorField&,
            const fvPatch&,
            const vectorField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new surfaceNormalFixedValueFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceNormalFixedValueFvPatchVectorField
        (
            const surfaceNormalFixedValueFvPatchVectorField&,
            const vectorField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const vectorField& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new surfaceNormalFixedValueFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }



    // Member functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
