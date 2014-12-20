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
    processorLduCoupledInterface

Description
    This class contains the generic elements of a processor interface

SourceFiles
    processorLduCoupledInterface.C

\*---------------------------------------------------------------------------*/

#ifndef processorLduCoupledInterface_H
#define processorLduCoupledInterface_H

#include "lduCoupledInterface.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorLduCoupledInterface Declaration
\*---------------------------------------------------------------------------*/

class processorLduCoupledInterface
:
    virtual public lduCoupledInterface
{

protected:

    // Protected member functions

        //- Raw field send function
        template<class Type>
        void send
        (
            const tmp<Field<Type> >&,
            const bool bufferdTransfer
        ) const;

        //- Raw field receive function returning field
        template<class Type>
        tmp<Field<Type> > receive(const label size) const;


        //- Raw field send function with data compression
        template<class Type>
        void compressedSend
        (
            const tmp<Field<Type> >&,
            const bool bufferdTransfer
        ) const;

        //- Raw field receive function with data compression returning field
        template<class Type>
        tmp<Field<Type> > compressedReceive(const label size) const;


public:

    //- Runtime type information
    TypeName("processorLduCoupledInterface");

    // Constructors

    // Destructor

        virtual ~processorLduCoupledInterface();


    // Member Functions

        // Access

            //- Return processor number
            virtual int myProcNo() const = 0;

            //- Return neigbour processor number
            virtual int neighbProcNo() const = 0;

            //- Is the transform required
            virtual bool doTransform() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;

        //- Transform given patch component field
        void transformProcCoupleField
        (
            scalarField& f,
            const direction cmpt
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorLduCoupledInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
