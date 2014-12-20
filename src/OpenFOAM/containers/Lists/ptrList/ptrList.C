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
    A ptrList<class T> is a 1D array of pointers to objects
    of T 'T', where the size of the array is known and used for
    subscript bounds checking, etc.
    The element operator [] returns a reference to the object
    rather than a pointer

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "ptrList.H"
#include "ptrListLoopM.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct with zero length

template<class T>
ptrList<T>::ptrList()
:
    ptrs_(),
    nextFree_(0)
{}


// Construct with length specified

template<class T>
ptrList<T>::ptrList(const label s)
:
    ptrs_(s, (T*)NULL),
    nextFree_(0)
{}


// Construct as copy

template<class T>
ptrList<T>::ptrList(const ptrList<T>& a)
:
    ptrs_(a.size()),
    nextFree_(a.nextFree_)
{
    forAll(*this, i)
    {
        ptrs_[i] = (a[i]).clone().ptr();
    }
}


// Construct as copy or re-use as specified.

template<class T>
ptrList<T>::ptrList(ptrList<T>& a, bool reUse)
:
    ptrs_(a.size()),
    nextFree_(a.nextFree_)
{
    if (reUse)
    {
        forAll(*this, i)
        {
            ptrs_[i] = a.ptrs_[i];
            a.ptrs_[i] = NULL;
        }
        a.nextFree_ = 0;
        a.setSize(0);
    }
    else
    {
        forAll(*this, i)
        {
            ptrs_[i] = (a[i]).clone().ptr();
        }
    }
}

// Construct as copy of SLPtrList<T>
template<class T>
ptrList<T>::ptrList(const SLPtrList<T>& sll)
:
    ptrs_(sll.size()),
    nextFree_(sll.size())
{
    if (sll.size())
    {
        label i = 0;
        for
        (
            typename SLPtrList<T>::const_iterator iter = sll.begin();
            iter != sll.end();
            ++iter
        )
        {
            ptrs_[i++] = (iter()).clone().ptr();
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy ptrlist elements

template<class T>
ptrList<T>::~ptrList()
{
    forAll(*this, i)
    {
        if (ptrs_[i])
        {
            delete ptrs_[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return an element pointer to be set.  The pointer is checked
// if already set and the index range checked.

template<class T>
typename ptrList<T>::Tptr& ptrList<T>::set(const label i)
{
    if (ptrs_[i])
    {
        FatalErrorIn("ptrList::set(const label)")
            << "pointer already set, cannot set to new element"
            << abort(FatalError);
    }

    return ptrs_[i];
}


// Hook an an element created by new onto the list

template<class T>
void ptrList<T>::hook(T* ptr)
{
    if (nextFree_ < size())
    {
        set(nextFree_) = ptr;
    }
    else
    {
        FatalErrorIn("ptrList::hook(T*)")
            << "ptrList full, cannot add new element"
            << abort(FatalError);
    }

    nextFree_++;
}


template<class T>
void ptrList<T>::setSize(const label newSize)
{
    label oldSize = size();

    if (newSize == 0)
    {
        clear();
    }
    else if (newSize < oldSize)
    {
        register label i;
        for (i=newSize; i<oldSize; i++)
        {
            if (ptrs_[i] != NULL)
            {
                FatalErrorIn("void ptrList<T>::setSize(const label newSize)")
                    << "sorry, cannot remove allocated ptrList elements"
                    << abort(FatalError);
            }
        }

        ptrs_.setSize(newSize);
    }
    else if (newSize > oldSize)
    {
        ptrs_.setSize(newSize);

        register label i;
        for (i=oldSize; i<newSize; i++)
        {
            ptrs_[i] = NULL;
        }
    }
}


template<class T>
void ptrList<T>::clear()
{
    forAll(*this, i)
    {
        if (ptrs_[i])
        {
            delete ptrs_[i];
        }
    }

    ptrs_.clear();

    // Reset the next free element to zero
    nextFree_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment

template<class T>
ptrList<T>& ptrList<T>::operator=(const ptrList<T>& a)
{
    if (this == &a)
    {
        FatalErrorIn("ptrList<T>::operator=(const ptrList<T>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (size() == 0)
    {
        setSize(a.size());
        nextFree_ = a.nextFree_;

        forAll(*this, i)
        {
            ptrs_[i] = (a[i]).clone().ptr();
        }
    }
    else if (a.size() == size())
    {
        forAll(*this, i)
        {
            (*this)[i] = a[i];
        }
    }
    else
    {
        FatalErrorIn("ptrList::operator=(const ptrList<T>&)")
            << "bad size: " << a.size()
            << abort(FatalError);
    }


    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ptrListIO.C"

// ************************************************************************* //
