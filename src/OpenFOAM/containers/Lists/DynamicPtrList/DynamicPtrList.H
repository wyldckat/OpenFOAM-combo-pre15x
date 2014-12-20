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
    DynamicPtrList

Description
    A dynamic list is a 1-D vector of objects of type T which resizes
    itself as necessary to accept the new objects.  Internal storage
    is a compact array and the list can be shrunk to compact storage.
    The increment of list size is controlled by two static members,
    which should be reset before the list is created.  This allows the
    list storage to either increate by the given increment or a given
    multiplier.
    This version is created to hold pointers to objects.

SourceFiles
    DynamicPtrListI.H
    DynamicPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicPtrList_H
#define DynamicPtrList_H

#include "ptrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T> class DynamicPtrList;

template<class T> Ostream& operator<<(Ostream&, const DynamicPtrList<T>&);

/*---------------------------------------------------------------------------*\
                         Class DynamicPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class DynamicPtrList
:
    public ptrList<T>
{
    // Static data members

        //- Static data to control resizing of the list
        static label sizeIncrement_;
        static label sizeMultiplier_;


public:

    // Constructors

        //- Construct null
        inline DynamicPtrList();

        //- Construct given size
        explicit inline DynamicPtrList(const label);


    // Member Functions

        // Access

            //- Size of the active part of the list.
            //  Direct over-ride of list size member function
            inline label size() const;


        // Edit

            //- Clear the list, i.e. set size to zero.
            void clear();

            //- Set the list size increment, returning the previous value
            inline static label setIncrement(const label);

            //- Set the list multiplier, returning the previous value
            inline static label setMultiplier(const label);

            //- Shrink the ptrList<T> to the number of elements used
            inline DynamicPtrList<T>& shrink();


    // Member Operators

        //- Append an element at the end of the list
        inline void append(const T& e);

        //- Return non-const access to an element,
        //  resizing the list if necessary
        inline T& operator()(const label);

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // Ostream operator

        // Write DynamicPtrList to Ostream.
        friend Ostream& operator<< <T>(Ostream&, const DynamicPtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicPtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
