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
    SubList

Description
    SubList is a List obtained as a section of another List.  Thus it is
    itself unallocated so that no storage is allocated or de-allocated
    during it's use.  To achieve this behaviour, SubList is derived from
    UList rather than List.

SourceFiles
    SubListI.H

\*---------------------------------------------------------------------------*/

#ifndef SubList_H
#define SubList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SubList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SubList
:
    public UList<T>
{
    // Private member functions

        //- Dissalow non-const return subscript-checked element of SubList.
        T& operator[](const label);


public:

    // Constructors

        //- Construct from UList and SubList size
        inline SubList
        (
            const UList<T>& list,
            const label subSize
        );

        //- Construct from UList, start and end indices
        inline SubList
        (
            const UList<T>& list,
            const label subSize,
            const label startIndex
        );


    // Member functions

        //- Return a null SubList
        static inline SubList<T>& null();


    // Member operators

        //- Return const subscript-checked element of SubList.
        const T& operator[](const label l) const
        {
            return UList<T>::operator[](l);
        }

        //- Allow cast to a const List<T>&
        inline operator const List<T>&() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
