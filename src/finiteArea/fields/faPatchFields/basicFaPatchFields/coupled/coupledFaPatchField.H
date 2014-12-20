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
    coupledFaPatchField

Description

SourceFiles
    coupledFaPatchField.C


\*---------------------------------------------------------------------------*/

#ifndef coupledFaPatchField_H
#define coupledFaPatchField_H

#include "faPatchField.H"
#include "coupledFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupledFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledFaPatchField
:
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("coupled");


    // Constructors

        //- Construct from patch and internal field
        coupledFaPatchField
        (
            const faPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        coupledFaPatchField
        (
            const faPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        coupledFaPatchField
        (
            const faPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping the given coupledFaPatchField onto a new patch
        coupledFaPatchField
        (
            const coupledFaPatchField<Type>&,
            const faPatch&,
            const Field<Type>&,
            const faPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        coupledFaPatchField
        (
            const coupledFaPatchField<Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone(const Field<Type>&) const = 0;


    // Member functions

        // Access

            //- Return true if this patch field is derived from
            //  coupledFaPatchField<Type>.
            virtual bool coupled() const
            {
                return true;
            }

            //- Initialise neighbour field transfer of given internal field
            virtual void initPatchNeighbourField(const Field<Type>&) const
            {}

            //- Return neighbour field given internal field
            virtual tmp<Field<Type> > patchNeighbourField
            (
                const Field<Type>&
            ) const = 0;


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad() const;

            //- Initialise neighbour field transfer of internal field
            virtual void initPatchNeighbourField() const;

            //- Return neighbour field of the internal field
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const bool);

            //- Evaluate the patch field
            virtual void evaluate();

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


        // Coupled interface functionality

            //- Return neighbour colouring
            virtual tmp<labelField> nbrColour
            (
                const labelField&
            ) const = 0;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction
            ) const = 0;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "coupledFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
