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

\*---------------------------------------------------------------------------*/

#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef USE_RANDOM

#   include <climits>

#   if INT_MAX    != 2147483647
#       error "INT_MAX    != 2147483647"
#       error "The random number generator random() may not work!"
#   endif

#else


#   ifdef cray
extern "C" {

     double drand48 (void);

     double erand48 (unsigned short xsubi[3]);

     long lrand48 (void);

     long nrand48 (unsigned short xsubi[3]);

     long mrand48 (void);

     long jrand48 (unsigned short xsubi[3]);

     void srand48 (long seedval);

     unsigned short *seed48 (unsigned short seed16v[3]);

     void lcong48 (unsigned short param[7]);
}
#   else

#       include <cstdlib>

#   endif

#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// construct given seed
Random::Random(const label& seed)
{
    if (seed > 1)
    {
        Seed = seed;
    }
    else
    {
        Seed = 1;
    }

#   ifdef USE_RANDOM
        srandom((unsigned int)Seed);
#   else
        srand48(Seed);
#   endif

}


int Random::bit()
{
#   ifdef USE_RANDOM
    if (random() > INT_MAX/2)
#   else
    if (mrand48() > 0)
#   endif
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


scalar Random::scalar01()
{
#   ifdef USE_RANDOM
        return (scalar)random()/INT_MAX;
#   else
        return drand48();
#   endif
}


vector Random::vector01()
{
    vector rndVec;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = scalar01();
    }

    return rndVec;
}


tensor Random::tensor01()
{
    tensor rndTen;
    for (direction cmpt=0; cmpt<tensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


vector Random::position(const vector& start, const vector& end)
{
    vector rndVec(start);

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) +=
            scalar01()*(end.component(cmpt) - start.component(cmpt));
    }

    return rndVec;
}


void Random::randomise(scalar& s)
{
     s = scalar01();
}


void Random::randomise(vector& v)
{
    v = vector01();
}


void Random::randomise(tensor& t)
{
    t = tensor01();
}


// return a normal Gaussian randon number
// with zero mean and unity variance N(0, 1)

scalar Random::GaussNormal()
{
    static int iset = 0;
    static scalar gset;
    scalar fac, rsq, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0*scalar01() - 1.0;
            v2 = 2.0*scalar01() - 1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;

        return v2*fac;
    }
    else
    {
        iset = 0;

        return gset;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
