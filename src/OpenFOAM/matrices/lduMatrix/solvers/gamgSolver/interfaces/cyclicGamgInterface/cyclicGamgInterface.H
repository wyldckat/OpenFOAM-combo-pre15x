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
    cyclicGamgInterface

Description
    GAMG agglomerated cyclic interface.

SourceFiles
    cyclicGamgInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGamgInterface_H
#define cyclicGamgInterface_H

#include "gamgInterface.H"
#include "cyclicLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGamgInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicGamgInterface
:
    public gamgInterface,
    virtual public cyclicLduInterface
{
    // Private data

        //- Reference tor the cyclicLduInterface from which this is
        //  agglomerated
        const cyclicLduInterface& fineCyclicInterface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicGamgInterface(const cyclicGamgInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicGamgInterface&);


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicGamgInterface
        (
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    // Destructor

        virtual ~cyclicGamgInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const unallocLabelList& iF
            ) const;


        //- Cyclic interface functions

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineCyclicInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return fineCyclicInterface_.reverseT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
