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
    processorFvPatchField

Description

SourceFiles
    processorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorFvPatchField_H
#define processorFvPatchField_H

#include "coupledFvPatchField.H"
#include "processorFvPatch.H"
#include "processorLduCoupledInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorFvPatchField
:
    virtual public processorLduCoupledInterface,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const processorFvPatch& procPatch_;


    // Private member functions

        //- Raw field send function
        template<class Type2>
        void send
        (
            const tmp<Field<Type2> >&,
            const bool bufferdTransfer
        ) const;

        //- Raw field receive function returning field
        template<class Type2>
        tmp<Field<Type2> > receive() const;


        //- Raw field send function with data compression
        template<class Type2>
        void compressedSend
        (
            const tmp<Field<Type2> >&,
            const bool bufferdTransfer
        ) const;

        //- Raw field receive function with data compression returning field
        template<class Type2>
        tmp<Field<Type2> > compressedReceive() const;


public:

    //- Runtime type information
    TypeName(processorFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorFvPatchField
        (
            const fvPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        processorFvPatchField
        (
            const fvPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorFvPatchField
        (
            const fvPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given processorFvPatchField onto a new patch
        processorFvPatchField
        (
            const processorFvPatchField<Type>&,
            const fvPatch&,
            const Field<Type>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        processorFvPatchField(const processorFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new processorFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        processorFvPatchField
        (
            const processorFvPatchField<Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone(const Field<Type>& iF) const
        {
            return tmp<fvPatchField<Type> >
            (
                new processorFvPatchField<Type>(*this, iF)
            );
        }


    // Destructor

        ~processorFvPatchField();


    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const
            {
                if (Pstream::parRun())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            //- Return neighbour field given internal field
            tmp<Field<Type> > patchNeighbourField() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const bool bufferdTransfer);

            //- Evaluate the patch field
            virtual void evaluate();

            //- Initialise neighbour colouring transfer
            virtual void initNbrColour
            (
                const labelField&,
                const bool bufferdTransfer
            ) const;

            //- Return neighbour colouring
            virtual tmp<labelField> nbrColour(const labelField&) const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const bool bufferdTransfer
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt
            ) const;

        //- Processor coupled interface functions

            //- Return processor number
            virtual int myProcNo() const
            {
                return procPatch_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return procPatch_.neighbProcNo();
            }

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return !(procPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procPatch_.forwardT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
