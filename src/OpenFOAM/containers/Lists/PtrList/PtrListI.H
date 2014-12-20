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

Description
    A PtrList<class T> is a 1D array of pointers to objects
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
inline label PtrList<T>::size() const
{
    return ptrs_.size();
}


template<class T>
inline void PtrList<T>::hook(const autoPtr<T>& aptr)
{
    hook(const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline void PtrList<T>::hook(const tmp<T>& t)
{
    hook(const_cast<tmp<T>&>(t).ptr());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Return element const reference

template<class T>
const T& PtrList<T>::operator[](const label i) const
{
    if (!ptrs_[i])
    {
        FatalErrorIn("PtrList::operator[] const")
            << "hanging pointer, cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


// Return element reference

template<class T>
T& PtrList<T>::operator[](const label i)
{
    if (!ptrs_[i])
    {
        FatalErrorIn("PtrList::operator[]")
            << "hanging pointer, cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
const T* PtrList<T>::operator()(const label i) const
{
    return ptrs_[i];
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline PtrList<T>::iterator::iterator(T** ptr)
:
    ptr_(ptr)
{}

template<class T>
inline bool PtrList<T>::iterator::operator==(const iterator& iter) const
{
    return ptr_ == iter.ptr_;
}

template<class T>
inline bool PtrList<T>::iterator::operator!=(const iterator& iter) const
{
    return ptr_ != iter.ptr_;
}

template<class T>
inline T& PtrList<T>::iterator::operator*()
{
    return **ptr_;
}

template<class T>
inline T& PtrList<T>::iterator::operator()()
{
    return operator*();
}

template<class T>
inline typename PtrList<T>::iterator
PtrList<T>::iterator::operator++()
{
    ++ptr_;
    return *this;
}

template<class T>
inline typename PtrList<T>::iterator
PtrList<T>::iterator::operator++(int)
{
    iterator tmp = *this;
    ++ptr_;
    return tmp;
}

template<class T>
inline typename PtrList<T>::iterator
PtrList<T>::iterator::operator--()
{
    --ptr_;
    return *this;
}

template<class T>
inline typename PtrList<T>::iterator
PtrList<T>::iterator::operator--(int)
{
    iterator tmp = *this;
    --ptr_;
    return tmp;
}

template<class T>
inline typename PtrList<T>::iterator
PtrList<T>::iterator::operator+=(label n)
{
    ptr_ += n;
    return *this;
}

template<class T>
inline typename PtrList<T>::iterator
operator+(const typename PtrList<T>::iterator& iter, label n)
{
    typename PtrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename PtrList<T>::iterator
operator+(label n, const typename PtrList<T>::iterator& iter)
{
    typename PtrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename PtrList<T>::iterator
PtrList<T>::iterator::operator-=(label n)
{
    ptr_ -= n;
    return *this;
}

template<class T>
inline typename PtrList<T>::iterator
operator-(const typename PtrList<T>::iterator& iter, label n)
{
    typename PtrList<T>::iterator tmp = iter;
    return tmp -= n;
}

template<class T>
inline label operator-
(
    const typename PtrList<T>::iterator& iter1,
    const typename PtrList<T>::iterator& iter2
)
{
    return (iter1.ptr_ - iter2.ptr_)/sizeof(T*);
}

template<class T>
inline T& PtrList<T>::iterator::operator[](label n)
{
    return *(*this + n);
}

template<class T>
inline bool PtrList<T>::iterator::operator<(const iterator& iter) const
{
    return ptr_ < iter.ptr_;
}

template<class T>
inline bool PtrList<T>::iterator::operator>(const iterator& iter) const
{
    return ptr_ > iter.ptr_;
}

template<class T>
inline bool PtrList<T>::iterator::operator<=(const iterator& iter) const
{
    return ptr_ <= iter.ptr_;
}

template<class T>
inline bool PtrList<T>::iterator::operator>=(const iterator& iter) const
{
    return ptr_ >= iter.ptr_;
}

template<class T>
inline typename PtrList<T>::iterator PtrList<T>::begin()
{
    return ptrs_.begin();
}

template<class T>
inline typename PtrList<T>::iterator PtrList<T>::end()
{
    return ptrs_.end();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
