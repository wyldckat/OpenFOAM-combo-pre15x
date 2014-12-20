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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "UList.H"
#include "ListLoopM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a null UList
template<class T>
UList<T>& UList<T>::null()
{
    UList<T>* nullPtr = (UList<T>*)NULL;
    return *nullPtr;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment of all entries to the given value
template<class T>
void UList<T>::operator=(const T& t)
{
    List_ACCESS(T, (*this), vp);
    List_FOR_ALL((*this), i)
        List_ELEM((*this), vp, i) = t;
    List_END_FOR_ALL
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T>
void UList<T>::swap(UList<T>& a)
{
    if (a.size_ != size_)
    {
        FatalErrorIn("UList<T>::swap(const UList<T>&)")
            << "ULists have different sizes"
            << abort(FatalError);
    }

    List_ACCESS(T, (*this), vp);
    List_ACCESS(T, a, ap);
    T tmp;
    List_FOR_ALL((*this), i)
        tmp = List_ELEM((*this), vp, i);
        List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
        List_ELEM(a, ap, i) = tmp;
    List_END_FOR_ALL
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Comparison for equality
template<class T>
bool UList<T>::operator==(const UList<T>& a) const
{
    if (size_ != a.size_)
    {
        return false;
    }

    bool equal = true;

    List_CONST_ACCESS(T, (*this), vp);
    List_CONST_ACCESS(T, (a), ap);

    List_FOR_ALL((*this), i)
        equal = equal && (List_ELEM((*this), vp, i) == List_ELEM((a), ap, i));
    List_END_FOR_ALL

    return equal;
}


// Comparison for inequality
template<class T>
bool UList<T>::operator!=(const UList<T>& a) const
{
    return !operator==(a);
}


// Compare ULists lexicographically
template<class T>
bool UList<T>::operator<(const UList<T>& a) const
{
    for
    (
        const_iterator vi = begin(), ai = a.begin();
        vi < end() && ai < a.end();
        vi++, ai++
    )
    {
        if (*vi < *ai)
        {
            return true;
        }
        else if (*vi > *ai)
        {
            return false;
        }
    }

    if (size_ < a.size_)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Compare ULists lexicographically
template<class T>
bool UList<T>::operator>(const UList<T>& a) const
{
    return a.operator<(*this);
}


//- Return true if !(a > b). Takes linear time.
template<class T>
bool UList<T>::operator<=(const UList<T>& a) const
{
    return !operator>(a);
}


//- Return true if !(a < b). Takes linear time.
template<class T>
bool UList<T>::operator>=(const UList<T>& a) const
{
    return !operator<(a);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "UListIO.C"

// ************************************************************************* //
