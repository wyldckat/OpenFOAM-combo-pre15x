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
    cyclicGamgInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    cyclicGamgInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGamgInterfaceField_H
#define cyclicGamgInterfaceField_H

#include "gamgInterfaceField.H"
#include "cyclicGamgInterface.H"
#include "cyclicLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGamgInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicGamgInterfaceField
:
    public gamgInterfaceField,
    virtual public cyclicLduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const cyclicGamgInterface& cyclicInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicGamgInterfaceField(const cyclicGamgInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicGamgInterfaceField&);
    

public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicGamgInterfaceField
        (
            const gamgInterface& gamgCp,
            const lduInterfaceField& fineInterfaceField
        );


    // Destructor

        virtual ~cyclicGamgInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicInterface_.size();
            }


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


        //- Cyclic interface functions

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicInterface_.reverseT();
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
