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
    DynamicList

Description
    A dynamic list is a 1-D vector of objects of type T which resizes
    itself as necessary to accept the new objects.  Internal storage
    is a compact array and the list can be shrunk to compact storage.
    The increment of list size is controlled by two template parameters,
    which allows the list storage to either increate by the given increment
    or a given multiplier.

SourceFiles
    DynamicListI.H
    DynamicList.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicList_H
#define DynamicList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
class DynamicList;

template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
Ostream& operator<<
(
    Ostream&,
    const DynamicList<T, SizeIncrement, SizeMultiplier>&
);


/*---------------------------------------------------------------------------*\
                           Class DynamicList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned SizeIncrement=0, unsigned SizeMultiplier=2>
class DynamicList
:
    public List<T>
{
    // Private data

        //- Number of next free element
        label nextFree_;


public:

    // Constructors

        //- Construct null
        inline DynamicList();

        //- Construct given size
        explicit inline DynamicList(const label);


    // Member Functions

        // Access

            //- Size of the active part of the list.
            //  Direct over-ride of list size member function
            inline label size() const;


        // Edit

            //- Reset size of List.
            void setSize(const label);

            //- Clear the list, i.e. set next free to zero.
            //  Allocated size does not change
            void clear();

            //- Shrink the List<T> to the number of elements used
            inline DynamicList<T, SizeIncrement, SizeMultiplier>& shrink();


    // Member Operators

        //- Append an element at the end of the list
        inline void append(const T& e);

        //- Return non-const access to an element,
        //  resizing the list if necessary
        inline T& operator()(const label);

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // Ostream operator

        // Write DynamicList to Ostream.
        friend Ostream& operator<< <T, SizeIncrement, SizeMultiplier>
        (
            Ostream&,
            const DynamicList<T, SizeIncrement, SizeMultiplier>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
