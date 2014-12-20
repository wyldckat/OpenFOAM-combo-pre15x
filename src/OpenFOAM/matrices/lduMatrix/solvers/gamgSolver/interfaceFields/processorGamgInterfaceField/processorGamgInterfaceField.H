/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    processorGamgInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    processorGamgInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef processorGamgInterfaceField_H
#define processorGamgInterfaceField_H

#include "gamgInterfaceField.H"
#include "processorGamgInterface.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorGamgInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class processorGamgInterfaceField
:
    public gamgInterfaceField,
    public processorLduInterfaceField
{
    // Private data

        //- Local reference cast into the processor interface
        const processorGamgInterface& procInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorGamgInterfaceField(const processorGamgInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const processorGamgInterfaceField&);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        processorGamgInterfaceField
        (
            const gamgInterface& gamgCp,
            const lduInterfaceField& fineInterface
        );


    // Destructor

        virtual ~processorGamgInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return procInterface_.size();
            }


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


        //- Processor interface functions

            //- Return processor number
            virtual int myProcNo() const
            {
                return procInterface_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return procInterface_.neighbProcNo();
            }

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procInterface_.forwardT();
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
