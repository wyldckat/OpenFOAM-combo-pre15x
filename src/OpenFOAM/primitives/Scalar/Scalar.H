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
    Scalar

Description
    Single floating point number

SourceFiles
    Scalar.C

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// template specialisation for pTraits<Scalar>
template<>
class pTraits<Scalar>
{
    Scalar p_;

public:

    //- Component type
    typedef Scalar cmptType;

    // Member constants

        enum
        {
            dim = 3,         // Dimensionality of space
            rank = 0,        // Rank od Scalar is 0
            nComponents = 1  // Number of components in Scalar is 1
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const Scalar zero;
        static const Scalar one;

    // Constructors

        //- Construct from Istream
        pTraits(Istream& is);

    // Member Functions

        operator Scalar() const
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return a string representation of a Scalar
word name(const Scalar s);


inline Scalar& setComponent(Scalar& s, const direction)
{
    return s;
}

inline Scalar component(const Scalar s, const direction)
{
    return s;
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2(y, x);
}

inline Scalar mag(const Scalar s)
{
    return ::fabs(s);
}

inline Scalar sign(const Scalar s)
{
    return (s >= 0)? 1: -1;
}

inline Scalar pos(const Scalar s)
{
    return (s >= 0)? 1: 0;
}

inline Scalar neg(const Scalar s)
{
    return (s < 0)? 1: 0;
}

inline bool equal(const Scalar& s1, const Scalar& s2)
{
    return mag(s1 - s2) <= ScalarVSMALL;
}

inline bool notEqual(const Scalar s1, const Scalar s2)
{
    return mag(s1 - s2) > ScalarVSMALL;
}

inline Scalar limit(const Scalar s1, const Scalar s2)
{
    return (mag(s1) < mag(s2))? s1: 0.0;
}

inline Scalar magSqr(const Scalar s)
{
    return s*s;
}

inline Scalar sqr(const Scalar s)
{
    return s*s;
}

inline Scalar pow3(const Scalar s)
{
    return s*s*s;
}

inline Scalar pow4(const Scalar s)
{
    return sqr(sqr(s));
}

inline Scalar cmptAv(const Scalar s)
{
    return s;
}

inline Scalar cmptMag(const Scalar s)
{
    return mag(s);
}

inline Scalar scale(const Scalar s, const Scalar d)
{
    return s*d;
}


#define transFunc(func)            \
inline Scalar func(const Scalar s) \
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

inline Scalar jn(const int n, const Scalar s)
{
    return ::jn(n, s);
}

transFunc(y0)
transFunc(y1)

inline Scalar yn(const int n, const Scalar s)
{
    return ::yn(n, s);
}

#undef transFunc


// Stabilisation around zero for division
inline Scalar stabilise(const Scalar s, const Scalar small)
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

Scalar readScalar(Istream& is);
Istream& operator>>(Istream&, Scalar&);
Ostream& operator<<(Ostream&, const Scalar);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
