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
    SortableList

Description

SourceFiles
    SortableList.C

\*---------------------------------------------------------------------------*/

#ifndef SortableList_H
#define SortableList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SortableList Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class SortableList
:
    public List<Type>
{
    // Private data

        //- Original indices
        labelList indices_;


public:

    // Public classes

        //- Less function class used by the sort function
        class less
        {
            const List<Type>& values_;

        public:

            less(const List<Type>& values)
            :
                values_(values)
            {}

            bool operator()(const label a, const label b)
            {
                return values_[a] < values_[b];
            }
        };


    // Constructors

        //- Construct from List, sorting the elements. Starts with indices set
        //  to index in argument
        SortableList(const List<Type>&);

        //- Construct from List, sorting the elements. Starts with indices 
        //  as given
        SortableList(const List<Type>&, const labelList&);

        //- Construct given size. Sort later on.
        SortableList(const label size);


    // Member Functions

        //- Return the list of sorted point indices
        const labelList& indices() const
        {
            return indices_;
        }

        //- Size the list (shrink or grow followed by sort)
        void setSize(const label);

        //- Sort the list (if changed after construction time)
        void sort();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SortableList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
