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
    fastLduTriLookup

Description
    The class implements fast off-diagonal lookup for the lduAddressing class.
    The problem here is to determine the index of the off-diagonal
    coefficient given the row and column index.  lduAddressing class
    provides the triIndex(const label a, const label b) which involves
    no setup cost but requires searching.

    During matrix assembly, it is necessary to use the triIndex(...) 
    member heavily and efficiency is crucial.  The fastLduTriLookup
    class implements the lookup without searching, but involves a
    setup and memory cost.  In order to do this, a subset of matrix
    addressing is be expended and a double direct lookup is used.

    Algorithm:
    The class allocates three sets of storage arrays:
    1) Diagonal lookup.  This is a diagonal array holding the index of
    the main addressing space where the row is expanded.  If the row
    is not expanded, the index is -1.  If a lookup is required for an
    inactive index, its row will be expanded in the next available
    slot.

    2) Main addressing space.  In this space, a number of arrays is
    expanded for direct lookup and is used for currently active rows.
    The size of main addressing space is configurable and determines
    the momery overhead.

    3) Reverse diagonal lookup.  When the main addressing space is
    re-used, the "old" diagonal index needs to know about it.  The
    reverse diagonal index gives the diagonal index for each row of
    the main addressing space.

    No clear-up of arrays is performed: the new row is expanded on top
    of the other and only the expanded values are used.

SourceFiles
    fastLduTriLookup.C

\*---------------------------------------------------------------------------*/

#ifndef fastLduTriLookup_H
#define fastLduTriLookup_H

#include "lduAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fastLduTriLookup Declaration
\*---------------------------------------------------------------------------*/

class fastLduTriLookup
{
    // Private data

        //- Reference to lduAddressing
        const lduAddressing& ldu_;

        //- Band (number of rows in main addressing space)
        const label band_;

        //- Diagonal addressing
        //  Row of main addressing space in which the row is expanded
        mutable labelList diagAddr_;

        //- Main addressing space
        //  Allocated as a single array and addressed using the offset
        mutable labelList mainAddr_;

        //- Reverse addressing
        //  Index of the diagonal for each row of the main addressing space
        mutable labelList reverseAddr_;

        //- Next free row of the main addressing space
        mutable label nextFree_;


        // Helper storage

            //- Low lookup index
            mutable label low_;

            //- High lookup index
            mutable label high_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fastLduTriLookup(const fastLduTriLookup&);

        //- Disallow default bitwise assignment
        void operator=(const fastLduTriLookup&);


        //- Expand a row into the next space and return index in the
        //  main addressing space
        void expandRow(const label d) const;


public:

    // Declare name of the class and its debug switch
    ClassName("fastLduTriLookup");

    // Constructors

        //- Construct from components
        fastLduTriLookup
        (
            const lduAddressing& ldu,
            const label band
        );


    // Destructor

        ~fastLduTriLookup();


    // Member Functions

        //- Return off-diagonal index given owner and neighbour label
        label triIndex(const label a, const label b) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
