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

#include "fvSchemes.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

int fvSchemes::debug(Foam::debug::debugSwitch("fvSchemes", false));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvSchemes::fvSchemes(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "fvSchemes",
            obr.time().system(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    ddtSchemes_(ITstream("ddtSchemes", tokenList())()),
    defaultDdtScheme_("default", tokenList()),
    d2dt2Schemes_(ITstream("d2dt2Schemes", tokenList())()),
    defaultD2dt2Scheme_("default", tokenList()),
    interpolationSchemes_(ITstream("interpolationSchemes", tokenList())()),
    defaultInterpolationScheme_("default", tokenList()),
    divSchemes_(ITstream("divSchemes", tokenList())()),
    defaultDivScheme_("default", tokenList()),
    gradSchemes_(ITstream("gradSchemes", tokenList())()),
    defaultGradScheme_("default", tokenList()),
    snGradSchemes_(ITstream("snGradSchemes", tokenList())()),
    defaultSnGradScheme_("default", tokenList()),
    laplacianSchemes_(ITstream("laplacianSchemes", tokenList())()),
    defaultLaplacianScheme_("default", tokenList()),
    fluxRequired_(ITstream("fluxRequired", tokenList())()),
    defaultFluxRequired_(false)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool fvSchemes::read()
{
    if (regIOobject::read())
    {
        const dictionary& dict = schemesDict();


        if (dict.found("ddtSchemes"))
        {
            ddtSchemes_ = dict.subDict("ddtSchemes");
        }
        else if (dict.found("timeScheme"))
        {
            // For backward compatibility.
            // The timeScheme will be deprecated with warning or removed in 2.4.

            word timeSchemeName(dict.lookup("timeScheme"));

            if (timeSchemeName == "EulerImplicit")
            {
                timeSchemeName = "Euler";
            }
            else if (timeSchemeName == "BackwardDifferencing")
            {
                timeSchemeName = "backward";
            }
            else if (timeSchemeName == "SteadyState")
            {
                timeSchemeName = "steadyState";
            }

            if (ddtSchemes_.found("default"))
            {
                ddtSchemes_.remove("default");
            }

            ddtSchemes_.add("default", timeSchemeName);
        }
        else
        {
            ddtSchemes_.add("default", "none");
        }

        if
        (
            ddtSchemes_.found("default")
         && word(ddtSchemes_.lookup("default")) != "none"
        )
        {
            defaultDdtScheme_ = ddtSchemes_.lookup("default");
        }


        if (dict.found("d2dt2Schemes"))
        {
            d2dt2Schemes_ = dict.subDict("d2dt2Schemes");
        }
        else if (dict.found("timeScheme"))
        {
            // For backward compatibility.
            // The timeScheme will be deprecated with warning or removed in 2.4.

            word timeSchemeName(dict.lookup("timeScheme"));

            if (timeSchemeName == "EulerImplicit")
            {
                timeSchemeName = "Euler";
            }
            else if (timeSchemeName == "SteadyState")
            {
                timeSchemeName = "steadyState";
            }

            if (d2dt2Schemes_.found("default"))
            {
                d2dt2Schemes_.remove("default");
            }

            d2dt2Schemes_.add("default", timeSchemeName);
        }
        else
        {
            d2dt2Schemes_.add("default", "none");
        }

        if
        (
            d2dt2Schemes_.found("default")
         && word(d2dt2Schemes_.lookup("default")) != "none"
        )
        {
            defaultD2dt2Scheme_ = d2dt2Schemes_.lookup("default");
        }


        if (dict.found("interpolationSchemes"))
        {
            interpolationSchemes_ = dict.subDict("interpolationSchemes");
        }
        else if (!interpolationSchemes_.found("default"))
        {
            interpolationSchemes_.add("default", "linear");
        }

        if
        (
            interpolationSchemes_.found("default")
         && word(interpolationSchemes_.lookup("default")) != "none"
        )
        {
            defaultInterpolationScheme_ =
                interpolationSchemes_.lookup("default");
        }


        divSchemes_ = dict.subDict("divSchemes");

        if
        (
            divSchemes_.found("default")
         && word(divSchemes_.lookup("default")) != "none"
        )
        {
            defaultDivScheme_ = divSchemes_.lookup("default");
        }


        gradSchemes_ = dict.subDict("gradSchemes");

        if
        (
            gradSchemes_.found("default")
         && word(gradSchemes_.lookup("default")) != "none"
        )
        {
            defaultGradScheme_ = gradSchemes_.lookup("default");
        }


        if (dict.found("snGradSchemes"))
        {
            snGradSchemes_ = dict.subDict("snGradSchemes");
        }
        else if (!snGradSchemes_.found("default"))
        {
            snGradSchemes_.add("default", "corrected");
        }

        if
        (
            snGradSchemes_.found("default")
         && word(snGradSchemes_.lookup("default")) != "none"
        )
        {
            defaultSnGradScheme_ = snGradSchemes_.lookup("default");
        }


        laplacianSchemes_ = dict.subDict("laplacianSchemes");

        if
        (
            laplacianSchemes_.found("default")
         && word(laplacianSchemes_.lookup("default")) != "none"
        )
        {
            defaultLaplacianScheme_ = laplacianSchemes_.lookup("default");
        }


        if (dict.found("fluxRequired"))
        {
            fluxRequired_ = dict.subDict("fluxRequired");

            if
            (
                fluxRequired_.found("default")
             && word(fluxRequired_.lookup("default")) != "none"
            )
            {
                defaultFluxRequired_ = Switch(fluxRequired_.lookup("default"));
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


const dictionary& fvSchemes::schemesDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


ITstream& fvSchemes::ddtScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup ddtScheme for " << name << endl;
    }

    if
    (
        ddtSchemes_.found(name)
     || !defaultDdtScheme_.size()
    )
    {
        return ddtSchemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultDdtScheme_)).rewind();
        return (ITstream&)defaultDdtScheme_;
    }
}


ITstream& fvSchemes::d2dt2Scheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup d2dt2Scheme for " << name << endl;
    }

    if
    (
        d2dt2Schemes_.found(name)
     || !defaultD2dt2Scheme_.size()
    )
    {
        return d2dt2Schemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultD2dt2Scheme_)).rewind();
        return (ITstream&)defaultD2dt2Scheme_;
    }
}


ITstream& fvSchemes::interpolationScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup interpolationScheme for " << name << endl;
    }

    if
    (
        interpolationSchemes_.found(name)
     || !defaultInterpolationScheme_.size()
    )
    {
        return interpolationSchemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultInterpolationScheme_)).rewind();
        return (ITstream&)defaultInterpolationScheme_;
    }
}


ITstream& fvSchemes::divScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup divScheme for " << name << endl;
    }

    if (divSchemes_.found(name) || !defaultDivScheme_.size())
    {
        return divSchemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultDivScheme_)).rewind();
        return (ITstream&)defaultDivScheme_;
    }
}


ITstream& fvSchemes::gradScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup gradScheme for " << name << endl;
    }

    if (gradSchemes_.found(name) || !defaultGradScheme_.size())
    {
        return gradSchemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultGradScheme_)).rewind();
        return (ITstream&)defaultGradScheme_;
    }
}


ITstream& fvSchemes::snGradScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup snGradScheme for " << name << endl;
    }

    if (snGradSchemes_.found(name) || !defaultSnGradScheme_.size())
    {
        return snGradSchemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultSnGradScheme_)).rewind();
        return (ITstream&)defaultSnGradScheme_;
    }
}


ITstream& fvSchemes::laplacianScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup laplacianScheme for " << name << endl;
    }

    if (laplacianSchemes_.found(name) || !defaultLaplacianScheme_.size())
    {
        return laplacianSchemes_.lookup(name);
    }
    else
    {
        ((ITstream&)(defaultLaplacianScheme_)).rewind();
        return (ITstream&)defaultLaplacianScheme_;
    }
}


bool fvSchemes::fluxRequired(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup fluxRequired for " << name << endl;
    }

    if (fluxRequired_.found(name))
    {
        return true;
    }
    else
    {
        return defaultFluxRequired_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
