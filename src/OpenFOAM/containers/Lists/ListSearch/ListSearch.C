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

Description

\*---------------------------------------------------------------------------*/

#include "ListSearch.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Find first occurence of given element and return index,
//  return -1 if not found
template<class List>
Foam::label Foam::findIndex(const List& l, typename List::const_reference t)
{
    label index = -1;

    forAll(l, i)
    {
        if (l[i] == t)
        {
            index = i;
            break;
        }
    }

    return index;
}


//- Find index of max element.
template<class List>
Foam::label Foam::findMax(const List& l)
{
    if (l.size() == 0)
    {
        return -1;
    }

    label index = 0;

    for (label i = 1; i < l.size(); i++)
    {
        if (l[i] > l[index])
        {
            index = i;
        }
    }

    return index;
}


//- Find index of min element.
template<class List>
Foam::label Foam::findMin(const List& l)
{
    if (l.size() == 0)
    {
        return -1;
    }

    label index = 0;

    for (label i = 1; i < l.size(); i++)
    {
        if (l[i] < l[index])
        {
            index = i;
        }
    }

    return index;
}


//- Find first occurence of given element in sorted list and return index,
//  return -1 if not found. Binary search.
template<class List>
Foam::label Foam::findSortedIndex
(
    const List& l,
    typename List::const_reference t
)
{
    if (l.size() == 0)
    {
        return -1;
    }

    label low = 0;
    label high = l.size() - 1;

    while (low <= high)
    {
        label mid = (low + high)/2;

        if (t < l[mid])
        {
            high = mid - 1;
        }
        else if (t > l[mid])
        {
            low = mid + 1;
        }
        else
        {
            return mid;
        }
    }

    return -1;
}


template<class List>
Foam::label Foam::findLower
(
    const List& l,
    typename List::const_reference t
)
{
    if (l.size() == 0)
    {
        return -1;
    }

    label low = 0;
    label high = l.size() - 1;

    while ((high - low) > 1)
    {
        label mid = (low + high)/2;

        if (l[mid] < t)
        {
            low = mid;
        }
        else
        {
            high = mid;
        }
    }

    if (l[high] < t)
    {
        return high;
    }
    else
    {
        if (l[low] < t)
        {
            return low;
        }
        else
        {
            return -1;
        }
    }
}


// ************************************************************************* //
