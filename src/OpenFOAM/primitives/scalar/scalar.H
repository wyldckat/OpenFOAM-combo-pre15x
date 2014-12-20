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

Type
    scalar

Description
    single floating point number identical to double

SourceFiles
    scalar.C

\*---------------------------------------------------------------------------*/

#ifndef scalar_H
#define scalar_H

#include "label.H"
#include "word.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define scalar as a double

namespace Foam
{
    typedef double scalar;

    // Largest and smallest scalar values allowed in certain parts of the code.
    // (15 is the number of significant figures in an
    //  IEEE double precision number.  See limits.h or float.h)
    static const scalar GREAT = 1.0e+15;
    static const scalar VGREAT = 1.0e+300;
    static const scalar SMALL = 1.0e-15;
    static const scalar VSMALL = 1.0e-300;
}

/*
// Define scalar as a float

namespace Foam
{
    typedef double scalar;

    // Largest and smallest scalar values allowed in certain parts of the code
    // (6 is the number of significant figures in an
    //  IEEE single precision number.  See limits.h or float.h)
    static const scalar GREAT = 1.0e+6;
    static const scalar VGREAT = 1.0e+37;
    static const scalar SMALL = 1.0e-6;
    static const scalar VSMALL = 1.0e-37;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pTraits.H"
#include "products.H"
#include "direction.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// template specialisation for pTraits<scalar>
template<>
class pTraits<scalar>
{
    scalar p_;

public:

    //- Component type
    typedef scalar cmptType;

    // Member constants

        enum
        {
            dim = 3,         // Dimensionality of space
            rank = 0,        // Rank od scalar is 0
            nComponents = 1  // Number of components in scalar is 1
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const scalar zero;
        static const scalar one;

    // Constructors

        //- Construct from Istream
        pTraits(Istream& is);

    // Member Functions

        operator scalar() const
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class typeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return a string representation of a scalar
word name(const scalar s);


inline scalar component(const scalar s, const direction)
{
    return s;
}


#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return ::pow(double(s), double(e));           \
}


MAXMINPOW(double, double, double)
MAXMINPOW(double, double, float)
MAXMINPOW(double, float, double)
MAXMINPOW(float, float, float)
MAXMINPOW(double, double, int)
MAXMINPOW(double, int, double)
MAXMINPOW(double, double, long)
MAXMINPOW(double, long, double)
MAXMINPOW(float, float, int)
MAXMINPOW(float, int, float)
MAXMINPOW(float, float, long)
MAXMINPOW(float, long, float)

#undef MAXMINPOW

inline scalar atan2(const scalar y, const scalar x)
{
    return ::atan2(y, x);
}

inline scalar mag(const scalar s)
{
    return ::fabs(s);
}

inline scalar sign(const scalar s)
{
    return (s >= 0)? 1: -1;
}

inline scalar pos(const scalar s)
{
    return (s >= 0)? 1: 0;
}

inline scalar neg(const scalar s)
{
    return (s < 0)? 1: 0;
}

template<class T>
inline bool equal(const T& s1, const T& s2)
{
    return s1 == s2;
}

inline bool equal(const scalar& s1, const scalar& s2)
{
    return mag(s1 - s2) <= VSMALL;
}

inline bool notEqual(const scalar s1, const scalar s2)
{
    return mag(s1 - s2) > VSMALL;
}

inline scalar limit(const scalar s1, const scalar s2)
{
    return (mag(s1) < mag(s2))? s1: 0.0;
}

inline scalar magSqr(const scalar s)
{
    return s*s;
}

inline scalar sqr(const scalar s)
{
    return s*s;
}

inline scalar pow3(const scalar s)
{
    return s*s*s;
}

inline scalar pow4(const scalar s)
{
    return sqr(sqr(s));
}

inline scalar cmptAv(const scalar s)
{
    return s;
}

inline scalar cmptMag(const scalar s)
{
    return mag(s);
}

inline scalar scale(const scalar s, const scalar d)
{
    return s*d;
}


#define transFunc(func)            \
inline scalar func(const scalar s) \
{                                  \
    return ::func(s);              \
}

// Standard C++ transcendental functions
transFunc(sqrt)
transFunc(exp)
transFunc(log)
transFunc(log10)
transFunc(sin)
transFunc(cos)
transFunc(tan)
transFunc(asin)
transFunc(acos)
transFunc(atan)
transFunc(sinh)
transFunc(cosh)
transFunc(tanh)
transFunc(asinh)
transFunc(acosh)
transFunc(atanh)

// Standard ANSI-C (but not in <cmath>) transcendental functions

transFunc(erf)
transFunc(erfc)
transFunc(lgamma)

transFunc(j0)
transFunc(j1)

inline scalar jn(const int n, const scalar s)
{
    return ::jn(n, s);
}

transFunc(y0)
transFunc(y1)

inline scalar yn(const int n, const scalar s)
{
    return ::yn(n, s);
}

#undef transFunc


// Stabilisation around zero for division
inline scalar stabilise(const scalar s, const scalar small)
{
    if (s >= 0)
    {
        return s + small;
    }
    else
    {
        return s - small;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

scalar readScalar(Istream& is);
Istream& operator>>(Istream&, scalar&);
Ostream& operator<<(Ostream&, const scalar);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
