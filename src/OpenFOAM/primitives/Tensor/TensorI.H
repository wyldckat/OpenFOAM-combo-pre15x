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
    Templated 3D Tensor derived from VectorSpace adding construction from
    9 components, element access using xx(), xy() etc. member functions and
    the iner-product (dot-product) and outer-product (tensor-product) operators.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template <class Cmpt>
inline Tensor<Cmpt>::Tensor()
{}


//- Construct given VectorSpace
template <class Cmpt>
inline Tensor<Cmpt>::Tensor(const VectorSpace<Tensor<Cmpt>, Cmpt, 9>& vs)
:
    VectorSpace<Tensor<Cmpt>, Cmpt, 9>(vs)
{}


//- Construct from components
template <class Cmpt>
inline Tensor<Cmpt>::Tensor
(
    const Cmpt txx, const Cmpt txy, const Cmpt txz,
    const Cmpt tyx, const Cmpt tyy, const Cmpt tyz,
    const Cmpt tzx, const Cmpt tzy, const Cmpt tzz
)
{
    this->v_[XX] = txx; this->v_[XY] = txy; this->v_[XZ] = txz;
    this->v_[YX] = tyx; this->v_[YY] = tyy; this->v_[YZ] = tyz;
    this->v_[ZX] = tzx; this->v_[ZY] = tzy; this->v_[ZZ] = tzz;
}


//- Construct from Istream
template <class Cmpt>
inline Tensor<Cmpt>::Tensor(Istream& is)
:
    VectorSpace<Tensor<Cmpt>, Cmpt, 9>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Cmpt>
inline const Vector<Cmpt>&  Tensor<Cmpt>::x() const
{
    return reinterpret_cast<const Vector<Cmpt>&>(this->v_[XX]);
}

template <class Cmpt>
inline const Vector<Cmpt>&  Tensor<Cmpt>::y() const
{
    return reinterpret_cast<const Vector<Cmpt>&>(this->v_[YX]);
}

template <class Cmpt>
inline const Vector<Cmpt>&  Tensor<Cmpt>::z() const
{
    return reinterpret_cast<const Vector<Cmpt>&>(this->v_[ZX]);
}


template <class Cmpt>
inline Vector<Cmpt>&  Tensor<Cmpt>::x()
{
    return reinterpret_cast<Vector<Cmpt>&>(this->v_[XX]);
}

template <class Cmpt>
inline Vector<Cmpt>&  Tensor<Cmpt>::y()
{
    return reinterpret_cast<Vector<Cmpt>&>(this->v_[YX]);
}

template <class Cmpt>
inline Vector<Cmpt>&  Tensor<Cmpt>::z()
{
    return reinterpret_cast<Vector<Cmpt>&>(this->v_[ZX]);
}


template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::xx() const
{
    return this->v_[XX];
}

template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::xy() const
{
    return this->v_[XY];
}

template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::xz() const
{
    return this->v_[XZ];
}


template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::yx() const
{
    return this->v_[YX];
}

template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::yy() const
{
    return this->v_[YY];
}

template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::yz() const
{
    return this->v_[YZ];
}


template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::zx() const
{
    return this->v_[ZX];
}

template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::zy() const
{
    return this->v_[ZY];
}

template <class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::zz() const
{
    return this->v_[ZZ];
}


template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::xx()
{
    return this->v_[XX];
}

template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::xy()
{
    return this->v_[XY];
}

template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::xz()
{
    return this->v_[XZ];
}


template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::yx()
{
    return this->v_[YX];
}

template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::yy()
{
    return this->v_[YY];
}

template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::yz()
{
    return this->v_[YZ];
}


template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::zx()
{
    return this->v_[ZX];
}

template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::zy()
{
    return this->v_[ZY];
}

template <class Cmpt>
inline Cmpt& Tensor<Cmpt>::zz()
{
    return this->v_[ZZ];
}


//- Return tensor transpose
template <class Cmpt>
inline Tensor<Cmpt> Tensor<Cmpt>::T() const
{
    return Tensor<Cmpt>
    (
        xx(), yx(), zx(),
        xy(), yy(), zy(),
        xz(), yz(), zz()
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Hodge Dual operator (tensor -> vector)
template <class Cmpt>
inline Vector<Cmpt> operator*(const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>(t.yz(), -t.xz(), t.xy());
}


//- Hodge Dual operator (vector -> tensor)
template <class Cmpt>
inline Tensor<Cmpt> operator*(const Vector<Cmpt>& v)
{
    return Tensor<Cmpt>
    (
             0, -v.z(),   v.y(),
         v.z(),      0,  -v.x(),
        -v.y(),  v.x(),       0
    );
}


//- Inner-product between two tensors
template <class Cmpt>
inline typename innerProduct<Tensor<Cmpt>, Tensor<Cmpt> >::type
operator&(const Tensor<Cmpt>& t1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*t2.xx() + t1.xy()*t2.yx() + t1.xz()*t2.zx(),
        t1.xx()*t2.xy() + t1.xy()*t2.yy() + t1.xz()*t2.zy(),
        t1.xx()*t2.xz() + t1.xy()*t2.yz() + t1.xz()*t2.zz(),

        t1.yx()*t2.xx() + t1.yy()*t2.yx() + t1.yz()*t2.zx(),
        t1.yx()*t2.xy() + t1.yy()*t2.yy() + t1.yz()*t2.zy(),
        t1.yx()*t2.xz() + t1.yy()*t2.yz() + t1.yz()*t2.zz(),

        t1.zx()*t2.xx() + t1.zy()*t2.yx() + t1.zz()*t2.zx(),
        t1.zx()*t2.xy() + t1.zy()*t2.yy() + t1.zz()*t2.zy(),
        t1.zx()*t2.xz() + t1.zy()*t2.yz() + t1.zz()*t2.zz()
    );
}

//- Inner-product between a tensor and a vector
template <class Cmpt>
inline typename innerProduct<Tensor<Cmpt>, Vector<Cmpt> >::type
operator&(const Tensor<Cmpt>& t, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        t.xx()*v.x() + t.xy()*v.y() + t.xz()*v.z(),
        t.yx()*v.x() + t.yy()*v.y() + t.yz()*v.z(),
        t.zx()*v.x() + t.zy()*v.y() + t.zz()*v.z()
    );
}

//- Inner-product between a vector and a tensor
template <class Cmpt>
inline typename innerProduct<Vector<Cmpt>, Tensor<Cmpt> >::type
operator&(const Vector<Cmpt>& v, const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>
    (
        v.x()*t.xx() + v.y()*t.yx() + v.z()*t.zx(),
        v.x()*t.xy() + v.y()*t.yy() + v.z()*t.zy(),
        v.x()*t.xz() + v.y()*t.yz() + v.z()*t.zz()
    );
}

//- Outer-product between two vectors
template <class Cmpt>
inline typename outerProduct<Vector<Cmpt>, Vector<Cmpt> >::type
operator*(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return Tensor<Cmpt>
    (
        v1.x()*v2.x(), v1.x()*v2.y(), v1.x()*v2.z(),
        v1.y()*v2.x(), v1.y()*v2.y(), v1.y()*v2.z(),
        v1.z()*v2.x(), v1.z()*v2.y(), v1.z()*v2.z()
    );
}


//- Return the diagonal of a tensor as a vector
template <class Cmpt>
inline Vector<Cmpt> diag(const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>(t.xx(), t.yy(), t.zz());
}

//- Return the trace of a tensor
template <class Cmpt>
inline scalar tr(const Tensor<Cmpt>& t)
{
    return t.xx() + t.yy() + t.zz();
}

//- Return the deviatoric part of a tensor
template <class Cmpt>
inline Tensor<Cmpt> dev(const Tensor<Cmpt>& t)
{
    return t - ((1.0/3.0)*tr(t))*Tensor<Cmpt>::I;
}

//- Return the deviatoric part of a tensor
template <class Cmpt>
inline Tensor<Cmpt> dev2(const Tensor<Cmpt>& t)
{
    return t - ((2.0/3.0)*tr(t))*Tensor<Cmpt>::I;
}

//- Return the symmetric part of a tensor
template <class Cmpt>
inline Tensor<Cmpt> symm(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        t.xx(), 0.5*(t.xy() + t.yx()), 0.5*(t.xz() + t.zx()),
        0.5*(t.yx() + t.xy()), t.yy(), 0.5*(t.yz() + t.zy()),
        0.5*(t.zx() + t.xz()), 0.5*(t.zy() + t.yz()), t.zz()
    );
}

//- Return the skew-symmetric part of a tensor
template <class Cmpt>
inline Tensor<Cmpt> skew(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        0.0, 0.5*(t.xy() - t.yx()), 0.5*(t.xz() - t.zx()),
        0.5*(t.yx() - t.xy()), 0.0, 0.5*(t.yz() - t.zy()),
        0.5*(t.zx() - t.xz()), 0.5*(t.zy() - t.yz()), 0.0
    );
}

//- Return the determinant of a tensor
template <class Cmpt>
inline scalar det(const Tensor<Cmpt>& t)
{
    return
    (
        t.xx()*t.yy()*t.zz() + t.xy()*t.yz()*t.zx()
      + t.xz()*t.yx()*t.zy() - t.xx()*t.yz()*t.zy()
      - t.xy()*t.yx()*t.zz() - t.xz()*t.yy()*t.zx()
    );
}

//- Return the cofactor tensor of a symmetric tensor
template <class Cmpt>
inline Tensor<Cmpt> cofactors(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        t.yy()*t.zz() - t.zy()*t.yz(),
        t.zx()*t.yz() - t.yx()*t.zz(),
        t.yx()*t.zy() - t.yy()*t.zx(),

        t.xz()*t.zy() - t.xy()*t.zz(),
        t.xx()*t.zz() - t.xz()*t.zx(),
        t.xy()*t.zx() - t.xx()*t.zy(),

        t.xy()*t.yz() - t.xz()*t.yy(),
        t.yx()*t.xz() - t.xx()*t.yz(),
        t.xx()*t.yy() - t.yx()*t.xy()
    );
}

//- Return the inverse of a symmetric tensor
template <class Cmpt>
inline Tensor<Cmpt> inv(const Tensor<Cmpt>& t)
{
    /*
    Tensor<Cmpt> t = s;

    if (mag(t.xx()) < SMALL)
    {
        t += Tensor<Cmpt>(1,0,0,0,0,0,0,0,0);
    }

    if (mag(t.yy()) < SMALL)
    {
        t += Tensor<Cmpt>(0,0,0,0,1,0,0,0,0);
    }

    if (mag(t.zz()) < SMALL)
    {
        t += Tensor<Cmpt>(0,0,0,0,0,0,0,0,1);
    }
    */
    return Tensor<Cmpt>
    (
        t.yy()*t.zz() - t.zy()*t.yz(),
        t.xz()*t.zy() - t.xy()*t.zz(),
        t.xy()*t.yz() - t.xz()*t.yy(),

        t.zx()*t.yz() - t.yx()*t.zz(),
        t.xx()*t.zz() - t.xz()*t.zx(),
        t.yx()*t.xz() - t.xx()*t.yz(),

        t.yx()*t.zy() - t.yy()*t.zx(),
        t.xy()*t.zx() - t.xx()*t.zy(),
        t.xx()*t.yy() - t.yx()*t.xy()
    )/det(t);
}


//- Return the 1st invariant of a tensor
template <class Cmpt>
inline scalar invariantI(const Tensor<Cmpt>& t)
{
    return tr(t);
}


//- Return the 2nd invariant of a tensor
template <class Cmpt>
inline scalar invariantII(const Tensor<Cmpt>& t)
{
    return
    (
        0.5*sqr(tr(t))
      - 0.5*
        (
           t.xx()*t.xx() + t.xy()*t.xy() + t.xz()*t.xz()
         + t.yx()*t.yx() + t.yy()*t.yy() + t.yz()*t.yz()
         + t.zx()*t.zx() + t.zy()*t.zy() + t.zz()*t.zz()
        )
    );
}


//- Return the 3rd invariant of a tensor
template <class Cmpt>
inline scalar invariantIII(const Tensor<Cmpt>& t)
{
    return det(t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
