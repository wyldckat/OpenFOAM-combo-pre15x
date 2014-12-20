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
    amgCoupledMatrixInterface

Description
     The coupled matrix multiplication for an lduMatrix is done by
     using the lduCoupledInterfaceField objects.  On the finest level,
     this role is performed by the patch fields.
     amgCoupledMatrixInterface performs the same function on coarse
     levels, as here the field is not available to perform the same
     function.


SourceFiles
    amgCoupledMatrixInterface.C

\*---------------------------------------------------------------------------*/

#ifndef amgCoupledMatrixInterface_H
#define amgCoupledMatrixInterface_H

#include "lduCoupledInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class amgCoupledMatrixInterface Declaration
\*---------------------------------------------------------------------------*/

class amgCoupledMatrixInterface
:
    public lduCoupledInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        amgCoupledMatrixInterface(const amgCoupledMatrixInterface&);

        //- Disallow default bitwise assignment
        void operator=(const amgCoupledMatrixInterface&);


public:

    // Constructors

        //- Construct from fine coupled interface
        amgCoupledMatrixInterface(const lduCoupledInterfaceField&);

    // Member Functions

        // Access

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt
            ) const
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
