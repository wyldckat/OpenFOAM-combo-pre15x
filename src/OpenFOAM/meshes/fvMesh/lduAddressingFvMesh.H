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
    lduAddressingFvMesh

Description

SourceFiles
    lduAddressingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef lduAddressingFvMesh_H
#define lduAddressingFvMesh_H

#include "lduAddressing.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class lduAddressingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class lduAddressingFvMesh
:
    public lduAddressing
{
    // Private data

        //- Lower as a subList of allOwner
        labelList::subList lowerAddr_;

        //- Upper as a reference to neighbour
        const labelList& upperAddr_;

        //- Patch addressing as a list of sublists
        PtrList<labelList::subList> patchAddr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduAddressingFvMesh(const lduAddressingFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduAddressingFvMesh&);


public:

    // Constructors

        //- Construct from components
        lduAddressingFvMesh(const fvMesh& m)
        :
            lduAddressing(m.nCells()),
            lowerAddr_
            (
                labelList::subList
                (
                    m.faceOwner(),
                    m.nInternalFaces()
                )
            ),
            upperAddr_(m.faceNeighbour()),
            patchAddr_(m.boundary().size())
        {
            forAll (m.boundary(), patchI)
            {
                patchAddr_.hook
                (
                    new labelList::subList(m.boundary()[patchI].faceCells())
                );
            }
        }


    // Destructor

        ~lduAddressingFvMesh()
        {}


    // Member Functions

        //- Return lower addressing (i.e. lower label = upper triangle)
        const unallocLabelList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
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
