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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "DLListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

DLListBase::iterator DLListBase::endIter
(
    const_cast<DLListBase&>(static_cast<const DLListBase&>(DLListBase())),
    reinterpret_cast<link*>(NULL)
);

DLListBase::const_iterator DLListBase::endConstIter
(
    static_cast<const DLListBase&>(DLListBase()),
    reinterpret_cast<const link*>(NULL)
);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// add to head of list
void DLListBase::insert(DLListBase::link* a)
{
    nElmts_++;

    if (!first_)
    {
        a->prev_ = 0;
        a->next_ = 0;
        first_ = last_ = a;
    }
    else
    {
        a->prev_ = 0;
        a->next_ = first_;
        first_->prev_ = a;
        first_ = a;
    }
}


// add to tail of list
void DLListBase::append(DLListBase::link* a)
{
    nElmts_++;

    if (!first_)
    {
        a->prev_ = 0;
        a->next_ = 0;
        first_ = last_ = a;
    }
    else
    {
        last_->next_ = a;
        a->prev_ = last_;
        a->next_ = 0;
        last_ = a;
    }
}


//- Swap this element with the one above unless it is at the top
bool DLListBase::swapUp(DLListBase::link* a)
{
    link* ap = a->prev_;

    if (ap)
    {
        if (ap == first_)
        {
            first_ = a;
        }

        if (a == last_)
        {
            last_ = ap;
        }

        if (a->next_)
        {
            a->next_->prev_ = ap;
        }

        if (ap->prev_)
        {
            ap->prev_->next_ = a;
        }

        a->prev_ = ap->prev_;
        ap->prev_ = a;

        ap->next_ = a->next_;
        a->next_ = ap;

        return true;
    }
    else
    {
        return false;
    }
}


//- Swap this element with the one below unless it is at the bottom
bool DLListBase::swapDown(DLListBase::link* a)
{
    link* an = a->next_;

    if (an)
    {
        if (a == first_)
        {
            first_ = an;
        }

        if (an == last_)
        {
            last_ = a;
        }

        if (a->prev_)
        {
            a->prev_->next_ = an;
        }

        if (an->next_)
        {
            an->next_->prev_ = a;
        }

        an->prev_ = a->prev_;
        a->prev_ = an;

        a->next_ = an->next_;
        an->next_ = a;

        return true;
    }
    else
    {
        return false;
    }
}


// remove and return head
DLListBase::link* DLListBase::removeHead()
{
    nElmts_--;

    if (!first_)
    {
        FatalErrorIn("void DLListBase::remove()")
            << "remove from empty list"
            << abort(FatalError);
    }

    DLListBase::link* f = first_;
    first_ = f->next_;

    if (!first_)
    {
        last_ = 0;
    }

    return f;
}


// remove and return element
DLListBase::link* DLListBase::remove(DLListBase::link* it)
{
    nElmts_--;

    link* ret = it;

    if (it == first_)
    {
        first_ = first_->next_;

        if (first_)
        {
            first_->prev_ = 0;
        }
    }
    else if (it == last_)
    {
        last_ = last_->prev_;

        if (last_)
        {
            last_->next_ = 0;
        }
    }
    else
    {
        it->next_->prev_ = it->prev_;
        it->prev_->next_ = it->next_;
    }

    return ret;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
