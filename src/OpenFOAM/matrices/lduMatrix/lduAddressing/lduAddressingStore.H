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
    lduAddressingStore

Description

\*---------------------------------------------------------------------------*/

#ifndef lduAddressingStore_H
#define lduAddressingStore_H

#include "lduAddressing.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduAddressingStore Declaration
\*---------------------------------------------------------------------------*/

class lduAddressingStore
:
    public lduAddressing
{
    // Private data

        //- Lower addressing
        labelList lowerAddr_;

        //- Upper addressing
        labelList upperAddr_;

        //- Patch to internal addressing
        labelListList patchAddr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduAddressingStore(const lduAddressingStore&);

        //- Disallow default bitwise assignment
        void operator=(const lduAddressingStore&);


public:

    // Constructors

        //- Construct from components
        lduAddressingStore
        (
            const label nEqns,
            const unallocLabelList& l,
            const unallocLabelList& u,
            const labelListList& pa
        )
        :
            lduAddressing(nEqns),
            lowerAddr_(l),
            upperAddr_(u),
            patchAddr_(pa)
        {}

    // Destructor

        ~lduAddressingStore()
        {}


    // Member Functions

        //- Return Lower addressing
        const unallocLabelList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return Upper addressing
        const unallocLabelList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return patch addressing
        const unallocLabelList& patchAddr(const label i) const
        {
            return patchAddr_[i];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
