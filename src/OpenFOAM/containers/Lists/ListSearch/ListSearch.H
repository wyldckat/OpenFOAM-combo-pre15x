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
    ListSearch

Description
    Various functions to search in Lists.

SourceFiles
    ListSearch.C

\*---------------------------------------------------------------------------*/

#ifndef ListSearch_H
#define ListSearch_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Find first occurence of given element and return index,
//  return -1 if not found. Linear search.
template<class List>
label findIndex(const List&, typename List::const_reference);


//- Find index of max element (and larger than given element).
//  return -1 if not found. Linear search.
template<class List>
label findMax(const List&);


//- Find index of min element (and less than given element).
//  return -1 if not found. Linear search.
template<class List>
label findMin(const List&);


//- Find first occurence of given element in sorted list and return index,
//  return -1 if not found. Binary search.
template<class List>
label findSortedIndex(const List&, typename List::const_reference);


//- Find last element < given value in sorted list and return index,
//  return -1 if not found. Binary search.
template<class List>
label findLower(const List&, typename List::const_reference);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ListSearch.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

