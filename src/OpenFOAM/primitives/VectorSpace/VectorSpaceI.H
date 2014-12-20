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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "products.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt>::VectorSpace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline const Cmpt& VectorSpace<Form, Cmpt, nCmpt>::component
(
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::component(direction) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt& VectorSpace<Form, Cmpt, nCmpt>::component
(
    const direction d
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn("VectorSpace<Form, Cmpt, nCmpt>::component(direction)")
            << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::component
(
    Cmpt& c,
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::component(Cmpt&, direction) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    c = v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::replace
(
    const direction d,
    const Cmpt& c
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::"
            "replace(direction, const Cmpt&) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    v_[d] = c;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline const Cmpt& VectorSpace<Form, Cmpt, nCmpt>::operator[]
(
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::operator[](direction d) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt& VectorSpace<Form, Cmpt, nCmpt>::operator[]
(
    const direction d
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn("VectorSpace<Form, Cmpt, nCmpt>::operator[](direction d)")
            << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator+=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    for (int i=0; i<nCmpt; i++)
    {
        v_[i] += vs.v_[i];
    }
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator-=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    for (int i=0; i<nCmpt; i++)
    {
        v_[i] -= vs.v_[i];
    }
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator*=
(
    const scalar s
)
{
    for (int i=0; i<nCmpt; i++)
    {
        v_[i] *= s;
    }
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator/=
(
    const scalar s
)
{
    for (int i=0; i<nCmpt; i++)
    {
        v_[i] /= s;
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline const Cmpt& component
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    const direction d
)
{
    return vs.component(d);
}


// Powers of a Form
// Equivalent to outer-products between the Form and itself

// Form^0 = 1.0
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 0>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>&,
    typename powProduct<Form, 0>::type
  = pTraits<typename powProduct<Form, 0>::type>::zero
)
{
    return 1.0;
}

// Form^1 = Form
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 1>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>& v,
    typename powProduct<Form, 1>::type
  = pTraits<typename powProduct<Form, 1>::type>::zero
)
{
    return static_cast<const Form&>(v);
}

// sqr(Form) = Form*Form
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 2>::type sqr
(
    const VectorSpace<Form, Cmpt, nCmpt>& v
)
{
    return static_cast<const Form&>(v)*static_cast<const Form&>(v);
}

// Form^2 = sqr(Form)
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 2>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>& v,
    typename powProduct<Form, 2>::type
  = pTraits<typename powProduct<Form, 2>::type>::zero
)
{
    return sqr(v);
}

/* Form^3 = Form*sqr(Form)
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 3>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>& v,
    typename powProduct<Form, 3>::type
  = pTraits<typename powProduct<Form, 3>::type>::zero
)
{
    return ((const Form&)v)*sqr(v);
}*/

/* Form^4 = sqr(sqr(Form))
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 4>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>& v,
    typename powProduct<Form, 4>::type
  = pTraits<typename powProduct<Form, 4>::type>::zero
)
{
    return sqr(sqr(v));
}*/


template<class Form, class Cmpt, int nCmpt>
inline scalar magSqr
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    scalar ms = magSqr(vs.v_[0]);
    for (int i=1; i<nCmpt; i++)
    {
        ms += magSqr(vs.v_[i]);
    }
    return ms;
}


template<class Form, class Cmpt, int nCmpt>
inline scalar mag
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    return ::sqrt(magSqr(vs));
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptMax
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt cMax = vs.v_[0];
    for (int i=1; i<nCmpt; i++)
    {
        cMax = max(cMax, vs.v_[i]);
    }
    return cMax;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptMin
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt cMin = vs.v_[0];
    for (int i=1; i<nCmpt; i++)
    {
        cMin = min(cMin, vs.v_[i]);
    }
    return cMin;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptSum
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt sum = vs.v_[0];
    for (int i=1; i<nCmpt; i++)
    {
        sum += vs.v_[i];
    }
    return sum;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptAv
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    return cmptSum(vs)/nCmpt;
}


template<class Form, class Cmpt, int nCmpt>
inline Form cmptMag
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = mag(vs.v_[i]);
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form max
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = max(vs1.v_[i], vs2.v_[i]);
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form min
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = min(vs1.v_[i], vs2.v_[i]);
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt> scale
(
    const VectorSpace<Form, Cmpt, nCmpt>& s,
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = scale(s.v_[i], vs.v_[i]);
    }
    return v;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline Form operator-
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = -vs.v_[i];
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator+
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = vs1.v_[i] + vs2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator-
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = vs1.v_[i] - vs2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator*
(
    scalar s,
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = s*vs.v_[i];
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator*
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    scalar s
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = vs.v_[i]*s;
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator/
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    scalar s
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = vs.v_[i]/s;
    }
    return v;
}

/*
template<class Form, class Cmpt, int nCmpt>
inline Form operator/
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = vs1.v_[i]/vs2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator/
(
    scalar s,
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    for (int i=0; i<nCmpt; i++)
    {
        v.v_[i] = s/vs.v_[i];
    }
    return v;
}
*/

template<class Form, class Cmpt, int nCmpt>
inline Cmpt operator&&
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Cmpt ddProd = vs1.v_[0]*vs2.v_[0];
    for (int i=1; i<nCmpt; i++)
    {
        ddProd += vs1.v_[i]*vs2.v_[i];
    }
    return ddProd;
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator==
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    bool eq = true;
    for (int i=0; i<nCmpt; i++)
    {
        if (!(eq &= (equal(vs1.v_[i], vs2.v_[i])))) break;
    }
    return eq;
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator!=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    return !(vs1 == vs2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
