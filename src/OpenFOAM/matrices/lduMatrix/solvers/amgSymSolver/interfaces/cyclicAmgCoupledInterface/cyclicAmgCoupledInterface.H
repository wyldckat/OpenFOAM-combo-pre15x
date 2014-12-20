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
    cyclicAmgCoupledInterface

Description

SourceFiles
    cyclicAmgCoupledInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAmgCoupledInterface_H
#define cyclicAmgCoupledInterface_H

#include "amgCoupledInterface.H"
#include "cyclicLduCoupledInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicAmgCoupledInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicAmgCoupledInterface
:
    public amgCoupledInterface,
    virtual public cyclicLduCoupledInterface
{
    // Private data

        //- Is the transform required
        bool doTransform_;

        // Face transformation tensor
        tensorField forwardT_;

        // Neighbour-cell transformation tensor
        tensorField reverseT_;

        //- Rank of component for transformation
        int rank_;

        //- Pointer to addressing
        mutable labelField* addrPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicAmgCoupledInterface(const cyclicAmgCoupledInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicAmgCoupledInterface&);

        //- Assemble a patch internal field
        template<class Type>
        tmp<Field<Type> > patchInternalField(const Field<Type>& iF) const;
    

public:

    //- Runtime type information
    TypeName("cyclic");

    // Constructors

        //- Construct from coupled interface
        cyclicAmgCoupledInterface
        (
            const lduCoupledInterface* fineInterfacePtr,
            const label index
        );

    // Member Functions

        // Access

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return true;
            }

            //- Return size
            label size() const;

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

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt
            ) const;

        //- Cyclic coupled interface functions

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

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return reverseT_;
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
