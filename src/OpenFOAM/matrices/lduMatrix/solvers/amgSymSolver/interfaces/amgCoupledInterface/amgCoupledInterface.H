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
    amgCoupledInterface

Description
     The coupled matrix multiplication for an lduMatrix is done by
     using the lduCoupledInterface objects.  On the finest level,
     this role is performed by the patch fields.
     amgCoupledInterface performs the same function on coarse
     levels, as here the field is not available to perform the same
     function.

SourceFiles
    amgCoupledInterface.C
    newAmgCoupledInterface.C

\*---------------------------------------------------------------------------*/

#ifndef amgCoupledInterface_H
#define amgCoupledInterface_H

#include "lduCoupledInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class amgCoupledInterface Declaration
\*---------------------------------------------------------------------------*/

class amgCoupledInterface
:
    virtual public lduCoupledInterface
{
    // Private data

        //- Index of the interface in the list
        const label index_;

        //- Dummy addressing array
        labelField addr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        amgCoupledInterface(const amgCoupledInterface&);

        //- Disallow default bitwise assignment
        void operator=(const amgCoupledInterface&);


public:

    //- Runtime type information
    TypeName("amgCoupledInterface");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            amgCoupledInterface,
            lduInterface,
            (const lduCoupledInterface* fineInterfacePtr, const label intIndex),
            (fineInterfacePtr, intIndex)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<amgCoupledInterface> New
        (
            const lduCoupledInterface* fineInterfacePtr,
            const label intIndex
        );


    // Constructors

        //- Construct from fine coupled interface
        amgCoupledInterface
        (
            const lduCoupledInterface*,
            const label intIndex
        )
        :
            index_(intIndex),
            addr_(0)
        {}


    // Member Functions

        // Access

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return interface index
            label index() const
            {
                return index_;
            }

            //- Return neighbour colouring
            virtual tmp<labelField> nbrColour
            (
                const labelField&
            ) const
            {
                return tmp<labelField>(new labelField(0));
            }

            //- Return interface coefficients given colouring and fine coeffs
            virtual tmp<scalarField> coeffs
            (
                const labelField&,
                const labelField&,
                const scalarField&
            ) const
            {
                // Dummy return for uncoupled patches
                return tmp<scalarField>(new scalarField(0));
            }

            //- Return addressing
            virtual const labelField& addressing() const
            {
                // Return dummy (zero-sized) addressing
                return addr_;
            }


        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceUpdate
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField& ,
                const direction
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction
            ) const
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
