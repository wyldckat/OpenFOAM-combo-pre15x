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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template <class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor()
{}


// Construct given VectorSpace
template <class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor
(
    const VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>& vs
)
:
    VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>(vs)
{}


// Construct given three Cmpts
template <class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor(const Cmpt& stii)
{
    this->v_[II] = stii;
}


// Construct from Istream
template <class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor(Istream& is)
:
    VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Cmpt>
inline const Cmpt&  SphericalTensor<Cmpt>::ii() const
{
    return this->v_[II];
}


template <class Cmpt>
inline Cmpt& SphericalTensor<Cmpt>::ii()
{
    return this->v_[II];
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template <class Cmpt>
inline Tensor<Cmpt>
operator+(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii() + t2.xx(), t2.xy(),            t2.xz(),
        t2.yx(),            st1.ii() + t2.yy(), t2.yz(),
        t2.zx(),            t2.zy(),            st1.ii() + t2.zz()
    );
}


template <class Cmpt>
inline Tensor<Cmpt>
operator+(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() + st2.ii(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() + st2.ii(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() + st2.ii()
    );
}


template <class Cmpt>
inline Tensor<Cmpt>
operator-(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii() - t2.xx(), -t2.xy(),            -t2.xz(),
       -t2.yx(),             st1.ii() - t2.yy(), -t2.yz(),
       -t2.zx(),            -t2.zy(),             st1.ii() - t2.zz()
    );
}


template <class Cmpt>
inline Tensor<Cmpt>
operator-(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() - st2.ii(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() - st2.ii(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() - st2.ii()
    );
}


//- Inner-product between two spherical tensors
template <class Cmpt>
inline SphericalTensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& st2)
{
    return SphericalTensor<Cmpt>(st1.ii()*st2.ii());
}


//- Inner-product between a spherical tensor and a tensor
template <class Cmpt>
inline Tensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii()*t2.xx(),
        st1.ii()*t2.xy(),
        st1.ii()*t2.xz(),

                          st1.ii()*t2.yx(),
                          st1.ii()*t2.yy(),
                          st1.ii()*t2.yz(),

                                            st1.ii()*t2.zx(),
                                            st1.ii()*t2.zy(),
                                            st1.ii()*t2.zz()
    );
}


//- Inner-product between a tensor and a spherical tensor
template <class Cmpt>
inline Tensor<Cmpt>
operator&(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*st2.ii(),
                          t1.xy()*st2.ii(),
                                            t1.xz()*st2.ii(),

        t1.yx()*st2.ii(),
                          t1.yy()*st2.ii(),
                                            t1.yz()*st2.ii(),

        t1.zx()*st2.ii(),
                          t1.zy()*st2.ii(),
                                            t1.zz()*st2.ii()
    );
}


//- Inner-product between a spherical tensor and a vector
template <class Cmpt>
inline Vector<Cmpt>
operator&(const SphericalTensor<Cmpt>& st, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        st.ii()*v.x(),
                       st.ii()*v.y(),
                                      st.ii()*v.z()
    );
}


//- Inner-product between a vector and a spherical tensor
template <class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const SphericalTensor<Cmpt>& st)
{
    return Vector<Cmpt>
    (
        v.x()*st.ii(),
                       v.y()*st.ii(),
                                      v.z()*st.ii()
    );
}


//- Division of a scalar by a sphericalTensor
template <class Cmpt>
inline SphericalTensor<Cmpt>
operator/(const scalar s, const SphericalTensor<Cmpt>& st)
{
    return SphericalTensor<Cmpt>(s/st.ii());
}


//- Return the trace of a spherical tensor
template <class Cmpt>
inline scalar tr(const SphericalTensor<Cmpt>& st)
{
    return 3*st.ii();
}


//- Return the determinant of a spherical tensor
template <class Cmpt>
inline scalar det(const SphericalTensor<Cmpt>& st)
{
    return st.ii()*st.ii()*st.ii();
}


//- Return the inverse of a symmetric tensor
template <class Cmpt>
inline SphericalTensor<Cmpt> inv(const SphericalTensor<Cmpt>& st)
{
    return SphericalTensor<Cmpt>(1.0/st.ii()); 
}


//- Return the deviatoric part of a tensor
template <class Cmpt>
inline Tensor<Cmpt> dev(const Tensor<Cmpt>& t)
{
    return t - SphericalTensor<Cmpt>::oneThirdI*tr(t);
}


//- Return the deviatoric part of a tensor
template <class Cmpt>
inline Tensor<Cmpt> dev2(const Tensor<Cmpt>& t)
{
    return t - SphericalTensor<Cmpt>::twoThirdsI*tr(t);
}


template<class Cmpt>
class typeOfSum<SphericalTensor<Cmpt>, Tensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};

template<class Cmpt>
class typeOfSum<Tensor<Cmpt>, SphericalTensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class outerProduct<SphericalTensor<Cmpt>, Cmpt>
{
public:

    typedef SphericalTensor<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SphericalTensor<Cmpt> >
{
public:

    typedef SphericalTensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, SphericalTensor<Cmpt> >
{
public:

    typedef SphericalTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, Tensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Tensor<Cmpt>, SphericalTensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, Vector<Cmpt> >
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SphericalTensor<Cmpt> >
{
public:

    typedef Vector<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
