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

Class
    messageStream

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "dictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int messageStream::level(debug::debugSwitch("level", 2));


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Construct from components
messageStream::messageStream
(
    const string& title,
    errorSeverity sev,
    const int maxErrors
)
:
    title_(title),
    severity_(sev),
    maxErrors_(maxErrors),
    errorCount_(0)
{}


//- Construct from dictionary
messageStream::messageStream(const dictionary& dict)
:
    title_(dict.lookup("title")),
    severity_(FATAL),
    maxErrors_(0),
    errorCount_(0)
{}


messageStream::operator OSstream&()
{
    if (level)
    {
        // Report the error
        if (severity_ == INFO && !Pstream::master())
        {
            return Snull;
        }
        else
        {
            if (title().size())
            {
                Sout<< title().c_str();
            }
        
            if (maxErrors_)
            {
                errorCount_++;
        
                if (errorCount_ >= maxErrors_)
                {
                    FatalErrorIn("messageStream::operator OSstream&()")
                        << "Too many errors"
                        << abort(FatalError);
                }
            }
        
            return Sout;
        }
    }

    return Snull;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global messageStream definitions

messageStream SeriousError 
(
    "--> FOAM Serious Error : ", messageStream::SERIOUS, 100
);

messageStream Warning("--> FOAM Warning : ", messageStream::WARNING);
messageStream Info("", messageStream::INFO);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
