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
    inkJetVelocityFvPatchVectorField

Description

SourceFiles
    inkJetVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef inkJetVelocityFvPatchVectorField_H
#define inkJetVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inkJetVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class inkJetVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        vector UBar_;
        vector UPrime_;
        scalar frequency_;


public:

    //- Runtime type information
    TypeName("inkJetVelocity");


    // Constructors

        //- Construct from patch and internal field
        inkJetVelocityFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&
        );

        //- Construct from patch, internal field and dictionary
        inkJetVelocityFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&,
            const dictionary&
        );

        //- Construct by mapping given inkJetVelocityFvPatchVectorField
        //  onto a new patch
        inkJetVelocityFvPatchVectorField
        (
            const inkJetVelocityFvPatchVectorField&,
            const fvPatch&,
            const vectorField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new inkJetVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        inkJetVelocityFvPatchVectorField
        (
            const inkJetVelocityFvPatchVectorField&,
            const vectorField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const vectorField& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new inkJetVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
