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
    lduCoupledInterface

Description
    A virtual base class containing the interface necessary to update
    an implicitly coupled interface in a lduMatrix.

\*---------------------------------------------------------------------------*/

#ifndef lduCoupledInterface_H
#define lduCoupledInterface_H

#include "labelField.H"
#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class lduMatrix;

/*---------------------------------------------------------------------------*\
                     Class lduCoupledInterface Declaration
\*---------------------------------------------------------------------------*/

class lduCoupledInterface
{
public:

    //- Runtime type information
    TypeName("lduCoupledInterface");


    // Constructors

        //- Construct null
        lduCoupledInterface()
        {}


    // Destructor

        virtual ~lduCoupledInterface();


    // Member Functions

        // Access

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Initialise neighbour colouring transfer
            virtual void initNbrColour
            (
                const labelField&,
                const bool
            ) const
            {}

            //- Return neighbour colouring
            virtual tmp<labelField> nbrColour
            (
                const labelField&
            ) const = 0;


        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const bool
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
            {
                notImplemented(type() + "::updateInterfaceMatrix()");
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
