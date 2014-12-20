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

Class
    error

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "OStringStream.H"
#include "fileName.H"
#include "dictionary.H"
#include "JobInfo.H"
#include "Pstream.H"
#include "OSspecific.H"

#if defined(__GNUC__)
#include "IStringStream.H"
#include <execinfo.h>
#include <demangle.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

error::error(const string& title)
:
    messageStream(title, messageStream::FATAL),
    functionName_("unknown"),
    sourceFileName_("unknown"),
    sourceFileLineNumber_(0),
    abort_(env("FOAM_ABORT")),
    throwExceptions_(false),
    messageStreamPtr_(new OStringStream())
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::error(const string& title) : cannot open error stream"
            << endl;
        ::exit(1);
    }
}


error::error(const dictionary& errDict)
:
    messageStream(errDict),
    functionName_(errDict.lookup("functionName")),
    sourceFileName_(errDict.lookup("sourceFileName")),
    sourceFileLineNumber_(readLabel(errDict.lookup("sourceFileLineNumber"))),
    abort_(env("FOAM_ABORT")),
    throwExceptions_(false),
    messageStreamPtr_(new OStringStream())
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::error(const dictionary& errDict) : "
               "cannot open error stream"
            << endl;
        ::exit(1);
    }
}


OSstream& error::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    messageStreamPtr_->rewind();
    functionName_ = functionName;
    sourceFileName_ = sourceFileName;
    sourceFileLineNumber_ = sourceFileLineNumber;

    return operator OSstream&();
}


OSstream& error::operator()
(
    const string& functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    return operator()
    (
        functionName.c_str(),
        sourceFileName,
        sourceFileLineNumber
    );
}


error::operator OSstream&()
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::operator OSstream&() : error stream has failed"
            << endl;
        printStack(Perr);
        ::abort();
    }

    return *messageStreamPtr_;
}


// Create and return a dictionary
error::operator dictionary() const
{
    dictionary errDict;

    string oneLineMessage(message());
    oneLineMessage.replaceAll('\n', ' ');

    errDict.add("type", word("Foam::error"));
    errDict.add("message", oneLineMessage);
    errDict.add("function", functionName());
    errDict.add("sourceFile", sourceFileName());
    errDict.add("sourceFileLineNumber", sourceFileLineNumber());

    return errDict;
}


string error::message() const
{
    return messageStreamPtr_->str();
}


void error::printStack(Ostream& os)
{
#if defined(__GNUC__)

    // Get raw stack symbols
    void *array[100];
    size_t size = backtrace(array, 100);
    char **strings = backtrace_symbols(array, size);

    // See if they contain function between () e.g. "(__libc_start_main+0xd0)"
    // and see if cplus_demangle can make sense of part before +
    for (size_t i = 0; i < size; i++)
    {
        string msg(strings[i]);

        string::size_type bracketPos = msg.find('(');

        if (bracketPos != string::size_type(string::npos))
        {
            string::size_type start = bracketPos+1;

            string::size_type plusPos = msg.find('+', start);

            if (plusPos != string::size_type(string::npos))
            {
                string cName(msg.substr(start, plusPos-start));

                char* cplusNamePtr = cplus_demangle
                (
                    cName.c_str(),
                    auto_demangling
                );

                if (cplusNamePtr)
                {
                    os<< cplusNamePtr << endl;
                    free(cplusNamePtr);
                }
                else
                {
                    os<< cName.c_str() << endl;
                }
            }
            else
            {
                string::size_type endBracketPos = msg.find(')', start);
                
                if (endBracketPos != string::size_type(string::npos))
                {
                    string fullName(msg.substr(start, endBracketPos-start));

                    os<< fullName.c_str() << endl;
                }
                else
                {
                    // Print raw message
                    os<< strings[i] << endl;
                }
            }
        }
        else
        {
            // Print raw message
            os<< strings[i] << endl;
        }
    }
    free(strings);

#endif
}


void error::exit(const int errNo)
{
    if (!throwExceptions_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        jobInfo.exit();
    }

    if (abort_)
    {
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        printStack(Perr);
        ::abort();
    }

    if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run exiting\n" << endl;
        Pstream::exit(errNo);
    }
    else
    {
        if (throwExceptions_)
        {
            throw *this;
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM exiting\n" << endl;
            ::exit(1);
        }
    }
}


void error::abort()
{
    if (!throwExceptions_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        jobInfo.abort();
    }

    if (abort_)
    {
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        printStack(Perr);
        ::abort();
    }

    if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run aborting\n" << endl;
        printStack(Perr);
        Pstream::abort();
    }
    else
    {
        if (throwExceptions_)
        {
            throw *this;
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM aborting\n" << endl;
            printStack(Perr);
            ::abort();
        }
    }
}


Ostream& operator<<(Ostream& os, const error& fErr)
{
    os  << endl << fErr.title().c_str()
        << fErr.message().c_str();

    if (error::level >= 2 && fErr.sourceFileLineNumber())
    {
        os  << endl << endl
            << "    From function " << fErr.functionName().c_str() << endl
            << "    in file " << fErr.sourceFileName().c_str()
            << " at line " << fErr.sourceFileLineNumber() << '.';
    }

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

error FatalError ("--> FOAM FATAL ERROR : ");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
