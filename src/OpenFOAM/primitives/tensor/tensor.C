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
    Tensor of scalars.

\*---------------------------------------------------------------------------*/

#include "tensor.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const tensor::typeName = "tensor";

template<>
const char* tensor::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const tensor tensor::zero
(
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
);

template<>
const tensor tensor::one
(
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
);

template<>
const tensor tensor::I
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vector eigenValues(const tensor& t)
{
    static const scalar small(1e-10);

    // return eigenvalues in ASCENDING ORDER OF ABSOLUTE VALUES !!!
    scalar a = -t.xx() - t.yy() - t.zz();

    scalar b = t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
             - t.xy()*t.yx() - t.xz()*t.zx() - t.yz()*t.zy();

    scalar c = - t.xx()*t.yy()*t.zz() - t.xy()*t.yz()*t.zx()
               - t.xz()*t.yx()*t.zy() + t.xz()*t.yy()*t.zx()
               + t.xy()*t.yx()*t.zz() + t.xx()*t.yz()*t.zy();

    scalar P = (a*a - 3.0*b)/9.0;
    scalar Q = (2*a*a*a - 9*a*b + 27*c)/54.0;

    scalar disc = Q*Q - P*P*P;

    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if
    (
        (
            mag(t.xy()) + mag(t.xz()) + mag(t.yx())
          + mag(t.yz()) + mag(t.zx()) + mag(t.zy())
        )
      < small
    )
    {
        // diagonal matrix
        i = t.xx();
        ii = t.yy();
        iii = t.zz();
    }
    else if (disc < -small)
    {
        // three different real roots
        scalar theta = acos(Q/::pow(P, 1.5));

        scalar sqrtP = sqrt(P);

        i = -2.0*sqrtP*cos(theta/3.0) - a/3.0;
        ii = -2.0*sqrtP*cos((theta + 2.0*mathematicalConstant::pi)/3.0) - a/3.0;
        iii = -2.0*sqrtP*cos((theta - 2.0*mathematicalConstant::pi)/3.0) - a/3.0;
    }
    else if (disc >= -small && disc <= small)
    {
        // zero discriminant
        if (mag(c) < small)
        {
            // there is a zero root; solve quadratic using algorithm from
            // Numerical Recepies, pp156
            scalar discQuad = sqr(b) - 4*a*c;

            if (discQuad >= -small)
            {
                scalar q = -0.5*(b + sign(b)*sqrt(mag(discQuad)));

                i = 0;
                ii = q/a;
                iii = c/q;
            }
            else
            {
                FatalErrorIn("eigenValues(const tensor&)")
                    << "complex eigenvalues in quadratic, discQuad: "
                    << discQuad << tab << "tensor: " << t
                    << abort(FatalError);
            }
        }
        else
        {
            if (Q >= -small && Q <= small)
            {
                // three equal real roots
                scalar root = -a/3;

                return vector(root, root, root);
            }
            else
            {
                // two equal real roots + 1 different root
                scalar r1f = -a/3 + sqrt(a*a - 3*b)/3;
                scalar r2f = -a -2*r1f;

                scalar r1s = -a/3 - sqrt(a*a - 3*b)/3;
                scalar r2s = -a -2*r1s;

                if (mag(2*r1f + r2f + a) < SMALL)
                {
                    i = r1f;
                    ii = r1f;
                    iii = r2f;
                }
                else
                {
                    i = r1s;
                    ii = r1s;
                    iii = r2s;
                }
            }
        }
    }
    else
    {
        // positive discriminant, complex roots
        FatalErrorIn("eigenValues(const tensor&)")
            << "complex eigenvalues detected, disc: " << disc << tab
            << "tensor: " << t
            << abort(FatalError);
    }


    // hard-coded shell sort
    if (mag(i) > mag(ii))
    {
        // swap
        scalar aux(i);
        i = ii;
        ii = aux;
    }

    if (mag(ii) > mag(iii))
    {
        // swap
        scalar aux(ii);
        ii = iii;
        iii = aux;
    }

    if (mag(i) > mag(ii))
    {
        // swap
        scalar aux(i);
        i = ii;
        ii = aux;
    }

    return vector(i, ii, iii);
}


vector eigenVector(const tensor& t, const scalar lambda)
{
    static const scalar small(1e-10);

    // construct the matrix.
    tensor lhs(t - lambda*I);

    vector y(vector::zero);     // start with a zero vector

    // If the first sub-determinant is non-zero, the eigen vector has a
    // in this direction component
    scalar firstSub = lhs.yy()*lhs.zz() - lhs.yz()*lhs.zy();
    scalar secondSub = -lhs.xx()*lhs.zz() + lhs.xz()*lhs.zx();
    scalar thirdSub = lhs.xx()*lhs.yy() - lhs.xy()*lhs.yx();

    if (firstSub > small || firstSub < -small)
    {
        y = vector
            (
                1.0,
                (-lhs.yx()*lhs.zz() - lhs.zx()*lhs.zy())/firstSub,
                (-lhs.yx()*lhs.yz() - lhs.zx()*lhs.yy())/firstSub
            );

        y /= mag(y);

        return y;
    }
    else if (secondSub > small || secondSub < -small)
    {
        y = vector
            (
                (-lhs.xy()*lhs.zz() - lhs.yz()*lhs.xz())/secondSub,
                1.0,
                (-lhs.xy()*lhs.zx() - lhs.yz()*lhs.xx())/secondSub
            );

        y /= mag(y);

        return y;
    }
    else if (thirdSub > small || thirdSub < -small)
    {
        y = vector
            (
                (-lhs.xz()*lhs.yy() - lhs.yz()*lhs.xz())/thirdSub,
                (-lhs.xz()*lhs.yx() - lhs.yz()*lhs.xx())/thirdSub,
                1.0
            );

        y /= mag(y);

        return y;
    }

    return y;
}


tensor eigenVectors(const tensor& t)
{
    // WARNING. Using the fact that the eigenvalues are in ascending
    // absolute order. Change at your peril.
    static const scalar small(1e-15);

    vector ev1(1, 0, 0);
    vector ev2(0, 1, 0);
    vector ev3(0, 0, 1);

    // check for singular eigenvector problem matrix
    if (mag(t.xy() + t.xz() + t.yz()) < small)
    {
        // hard-coded shell sort
        if (mag(t.xx()) > mag(t.yy()))
        {
            // swap ev1 and ev2
            vector aux(ev1);
            ev1 = ev2;
            ev2 = aux;
        }
        if (max(mag(t.xx()), mag(t.yy())) > mag(t.zz()))
        {
            // swap ev2 and ev3
            vector aux(ev2);
            ev2 = ev3;
            ev3 = aux;
        }
        if (min(mag(t.xx()), mag(t.yy())) > mag(t.zz()))
        {
            // swap ev1 and ev2
            vector aux(ev1);
            ev1 = ev2;
            ev2 = aux;
        }
    }
    else
    {
        // solve for eigen vectors

        vector evals(eigenValues(t));

        ev1 = eigenVector(t, evals.x());
        ev2 = eigenVector(t, evals.y());
        ev3 = eigenVector(t, evals.z());
    }

    return tensor
    (
        ev1.x(), ev1.y(), ev1.z(),
        ev2.x(), ev2.y(), ev2.z(),
        ev3.x(), ev3.y(), ev3.z()
    );
}


// Matrix inversion with singular value decomposition
tensor hinv(const tensor& t)
{
    static const scalar large(1e10);
    static const scalar small(1e-10);

    if (det(t) > small)
    {
        return inv(t);
    }
    else
    {
        vector eig = eigenValues(t);
        tensor eigVecs = eigenVectors(t);

        tensor zeroInv(tensor::zero);

        if (mag(eig.z()) > large*mag(eig.x()))
        {
            zeroInv += sqr(vector(eigVecs.xx(), eigVecs.xy(), eigVecs.xz()));
        }

        if (mag(eig.z()) > large*mag(eig.y()))
        {
            // singular direction 1
            zeroInv += sqr(vector(eigVecs.yx(), eigVecs.yy(), eigVecs.yz()));
        }

        return inv(t + zeroInv) - zeroInv;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
