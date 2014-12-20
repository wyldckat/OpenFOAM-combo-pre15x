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
    Signal handler for floating point exceptions, abort is called.

\*---------------------------------------------------------------------------*/

#include "JobInfo.H"
#include "OSspecific.H"

#if defined(linux) || defined(linuxAMD64) || defined(linuxIA64) && defined(__GNUC__)

#   ifndef __USE_GNU
#       define __USE_GNU
#   endif

#   include <fenv.h>

#elif defined(sgiN32) || defined(sgiN32Gcc)

#   include <sigfpe.h>

#endif

#include <cstdlib>
#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigfpe Declaration
\*---------------------------------------------------------------------------*/

class sigfpe
{

public:

#   if defined(linux) || defined(linuxAMD64) || defined(linuxIA64) && defined(__GNUC__)

    static void sigfpeHandler(int)
    {
        Info << "Caught arithmetic exception" << endl;

        if (!jobInfo.found("termination"))
        {
            jobInfo.add("termination", word("SIGFPE"));
        }

        ::abort();
    }

#   endif

    sigfpe()
    {
        if (env("FOAM_SIGFPE"))
        {
#           if defined(linux) || defined(linuxAMD64) || defined(linuxIA64) && defined(__GNUC__)

            feenableexcept
            (
                FE_DIVBYZERO
              | FE_INVALID
              | FE_OVERFLOW
            );

            signal(SIGFPE, sigfpeHandler);

#           elif defined(sgiN32) || defined(sgiN32Gcc)

            sigfpe_[_DIVZERO].abort=1;
            sigfpe_[_OVERFL].abort=1;
            sigfpe_[_INVALID].abort=1;

            sigfpe_[_DIVZERO].trace=1;
            sigfpe_[_OVERFL].trace=1;
            sigfpe_[_INVALID].trace=1;

            handle_sigfpes
            (
                _ON,
                _EN_DIVZERO
              | _EN_INVALID
              | _EN_OVERFL,
                0,
                _ABORT_ON_ERROR,
                NULL
            );

#           endif
        }
    }
};


static const sigfpe initSigfpe;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
