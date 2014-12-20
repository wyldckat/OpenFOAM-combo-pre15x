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
    supersonicFreeStreamFvPatchVectorField

Description
    Supersonic free stream condition.

    Supersonic outflow is vented according to ???

    Supersonic inflow is assumed to occur according to the Prandtl-Meyer
    expansion process.

    Subsonic outflow is zero-gradiented from inside the domain.

    N.B. This boundary condition is ill-posed if the free-stream flow is
         normal to the boundary.

SourceFiles
    supersonicFreeStreamFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef supersonicFreeStreamFvPatchVectorFields_H
#define supersonicFreeStreamFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class supersonicFreeStreamFvPatch Declaration
\*---------------------------------------------------------------------------*/

class supersonicFreeStreamFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Velocity of the free stream
        vector UInf_;

        //- Pressure of the free stream
        scalar pInf_;

        //- Temperature of the free stream
        scalar TInf_;


public:

    //- Runtime type information
    TypeName("supersonicFreeStream");


    // Constructors

        //- Construct from patch and internal field
        supersonicFreeStreamFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&
        );

        //- Construct from patch, internal field and dictionary
        supersonicFreeStreamFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&,
            const dictionary&
        );

        //- Construct by mapping given supersonicFreeStreamFvPatchVectorField
        //  onto a new patch
        supersonicFreeStreamFvPatchVectorField
        (
            const supersonicFreeStreamFvPatchVectorField&,
            const fvPatch&,
            const vectorField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new supersonicFreeStreamFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        supersonicFreeStreamFvPatchVectorField
        (
            const supersonicFreeStreamFvPatchVectorField&,
            const vectorField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const vectorField& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new supersonicFreeStreamFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the velocity at infinity
            const vector& UInf() const
            {
                return UInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            vector& UInf()
            {
                return UInf_;
            }

            //- Return the pressure at infinity
            scalar pInf() const
            {
                return pInf_;
            }

            //- Return reference to the pressure at infinity to allow adjustment
            scalar& pInf()
            {
                return pInf_;
            }

            //- Return the temperature at infinity
            scalar TInf() const
            {
                return TInf_;
            }

            //- Return reference to the temperature at infinity 
            //  to allow adjustment
            scalar& TInf()
            {
                return TInf_;
            }


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
