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
    rather than a pointer.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline label ptrList<T>::size() const
{
    return ptrs_.size();
}


template<class T>
inline void ptrList<T>::hook(const autoPtr<T>& aptr)
{
    hook(((autoPtr<T>&)aptr).ptr());
}


template<class T>
inline void ptrList<T>::hook(const tmp<T>& t)
{
    hook(((tmp<T>&)t).ptr());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Return element const reference

template<class T>
const T& ptrList<T>::operator[](const label i) const
{
    if (!ptrs_[i])
    {
        FatalErrorIn("ptrList::operator[] const")
            << "hanging pointer, cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


// Return element reference

template<class T>
T& ptrList<T>::operator[](const label i)
{
    if (!ptrs_[i])
    {
        FatalErrorIn("ptrList::operator[]")
            << "hanging pointer, cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
const T* ptrList<T>::operator()(const label i) const
{
    return ptrs_[i];
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline ptrList<T>::iterator::iterator(T** ptr)
:
    ptr_(ptr)
{}

template<class T>
inline bool ptrList<T>::iterator::operator==(const iterator& iter) const
{
    return ptr_ == iter.ptr_;
}

template<class T>
inline bool ptrList<T>::iterator::operator!=(const iterator& iter) const
{
    return ptr_ != iter.ptr_;
}

template<class T>
inline T& ptrList<T>::iterator::operator*()
{
    return **ptr_;
}

template<class T>
inline T& ptrList<T>::iterator::operator()()
{
    return operator*();
}

template<class T>
inline typename ptrList<T>::iterator
ptrList<T>::iterator::operator++()
{
    ++ptr_;
    return *this;
}

template<class T>
inline typename ptrList<T>::iterator
ptrList<T>::iterator::operator++(int)
{
    iterator tmp = *this;
    ++ptr_;
    return tmp;
}

template<class T>
inline typename ptrList<T>::iterator
ptrList<T>::iterator::operator--()
{
    --ptr_;
    return *this;
}

template<class T>
inline typename ptrList<T>::iterator
ptrList<T>::iterator::operator--(int)
{
    iterator tmp = *this;
    --ptr_;
    return tmp;
}

template<class T>
inline typename ptrList<T>::iterator
ptrList<T>::iterator::operator+=(label n)
{
    ptr_ += n;
    return *this;
}

template<class T>
inline typename ptrList<T>::iterator
operator+(const typename ptrList<T>::iterator& iter, label n)
{
    typename ptrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename ptrList<T>::iterator
operator+(label n, const typename ptrList<T>::iterator& iter)
{
    typename ptrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename ptrList<T>::iterator
ptrList<T>::iterator::operator-=(label n)
{
    ptr_ -= n;
    return *this;
}

template<class T>
inline typename ptrList<T>::iterator
operator-(const typename ptrList<T>::iterator& iter, label n)
{
    typename ptrList<T>::iterator tmp = iter;
    return tmp -= n;
}

template<class T>
inline label operator-
(
    const typename ptrList<T>::iterator& iter1,
    const typename ptrList<T>::iterator& iter2
)
{
    return (iter1.ptr_ - iter2.ptr_)/sizeof(T*);
}

template<class T>
inline T& ptrList<T>::iterator::operator[](label n)
{
    return *(*this + n);
}

template<class T>
inline bool ptrList<T>::iterator::operator<(const iterator& iter) const
{
    return ptr_ < iter.ptr_;
}

template<class T>
inline bool ptrList<T>::iterator::operator>(const iterator& iter) const
{
    return ptr_ > iter.ptr_;
}

template<class T>
inline bool ptrList<T>::iterator::operator<=(const iterator& iter) const
{
    return ptr_ <= iter.ptr_;
}

template<class T>
inline bool ptrList<T>::iterator::operator>=(const iterator& iter) const
{
    return ptr_ >= iter.ptr_;
}

template<class T>
inline typename ptrList<T>::iterator ptrList<T>::begin()
{
    return ptrs_.begin();
}

template<class T>
inline typename ptrList<T>::iterator ptrList<T>::end()
{
    return ptrs_.end();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
