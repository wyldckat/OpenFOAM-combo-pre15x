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
    processorAmgCoupledInterface

Description

SourceFiles
    processorAmgCoupledInterface.C

\*---------------------------------------------------------------------------*/

#ifndef processorAmgCoupledInterface_H
#define processorAmgCoupledInterface_H

#include "amgCoupledInterface.H"
#include "processorLduCoupledInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorAmgCoupledInterface Declaration
\*---------------------------------------------------------------------------*/

class processorAmgCoupledInterface
:
    public amgCoupledInterface,
    virtual public processorLduCoupledInterface
{
    // Private data

        //- Local processor number
        int myProcNo_;

        //- Neighbour processor number
        int neighbProcNo_;

        //- Is the transform required
        bool doTransform_;

        // Face transformation tensor
        tensorField forwardT_;

        //- Rank of component for transformation
        int rank_;

        //- Pointer to addressing
        mutable labelField* addrPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorAmgCoupledInterface(const processorAmgCoupledInterface&);

        //- Disallow default bitwise assignment
        void operator=(const processorAmgCoupledInterface&);

        //- Assemble a patch internal field
        template<class Type>
        tmp<Field<Type> > patchInternalField(const Field<Type>& iF) const;
    

public:

    //- Runtime type information
    TypeName("processor");

    // Constructors

        //- Construct from coupled interface
        processorAmgCoupledInterface
        (
            const lduCoupledInterface* fineInterfacePtr,
            const label index
        );


    // Destructor

        virtual ~processorAmgCoupledInterface();


    // Member Functions

        // Access

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return true;
            }

            //- Return size
            label size() const;

            //- Initialise neighbour colouring transfer
            virtual void initNbrColour(const labelField&, const bool) const;

            //- Return neighbour colouring
            virtual tmp<labelField> nbrColour
            (
                const labelField&
            ) const;

            //- Given colouring for both sides, return coarse level coefficients
            //  Because of the algorithmic complexity, this also
            //  calculates the addressing.  
            virtual tmp<scalarField> coeffs
            (
                const labelField& localColour,
                const labelField& nbrColour,
                const scalarField& fineCoeffs
            ) const;

            //- Return addressing
            virtual const labelField& addressing() const;


        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const bool bufferedTransfer
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt
            ) const;


        //- Processor coupled interface functions

            //- Return processor number
            virtual int myProcNo() const
            {
                return myProcNo_;
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return neighbProcNo_;
            }

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
