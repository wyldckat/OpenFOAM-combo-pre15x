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

#include "error.H"

#include "faSchemes.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faSchemes::faSchemes(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "faSchemes",
            obr.time().system(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    timeScheme_(EI),
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
    fluxRequired_(ITstream("fluxRequired", tokenList())())
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool faSchemes::read()
{
    if (regIOobject::read())
    {
        const dictionary& dict = schemesDict();

        word timeSchemeName(dict.lookup("timeScheme"));
        
        if (timeSchemeName == "SteadyState")
        {
            timeScheme_ = SS;
        }
        else if (timeSchemeName == "EulerImplicit")
        {
            timeScheme_ = EI;
        }
        else if (timeSchemeName == "BackwardDifferencing")
        {
            timeScheme_ = BD;
        }
        else if (timeSchemeName == "CrankNicholson")
        {
            timeScheme_ = CN;

            FatalIOErrorIn("faSchemes::read()", *this)
                << timeSchemeName << " is not currently supported"
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorIn("faSchemes::read()", *this)
                << "unknown timeScheme "
                << timeSchemeName << endl
                << "    Should be one of "
                   "SteadyState"
                   ", EulerImplicit"
                   " or BackwardDifferencing"
                << exit(FatalIOError);
        }

        if (dict.found("interpolationSchemes"))
        {
            interpolationSchemes_ = dict.subDict("interpolationSchemes");
        }
        else
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
        else
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
        }

        return true;
    }
    else
    {
        return false;
    }
}


const dictionary& faSchemes::schemesDict() const
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


faSchemes::timeSchemes faSchemes::timeScheme() const
{
    return timeScheme_;
}


ITstream& faSchemes::interpolationScheme(const word& name) const
{
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


ITstream& faSchemes::divScheme(const word& name) const
{
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


ITstream& faSchemes::gradScheme(const word& name) const
{
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


ITstream& faSchemes::snGradScheme(const word& name) const
{
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


ITstream& faSchemes::laplacianScheme(const word& name) const
{
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


bool faSchemes::fluxRequired(const word& name) const
{
    return fluxRequired_.found(name);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
