/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "error.H"
#include "OStringStream.H"
#include "fileName.H"
#include "dictionary.H"
#include "JobInfo.H"
#include "Pstream.H"
#include "OSspecific.H"

#if defined(__GNUC__)
#if defined(__linux__) && (defined(__i386__) || defined(__x86_64__))
#define HAS_DEMANGLING
#endif
#endif

#ifdef HAS_DEMANGLING
#include "IStringStream.H"
#include "IFstream.H"
#include "readHexLabel.H"
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


error::~error()
{
    delete messageStreamPtr_;
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



#if defined(HAS_DEMANGLING)

string doThePOpenDance(const string &cmd, label line=0)
{
    const int MAX=1000;

    FILE *cmdPipe = popen(cmd.c_str(),"r");
    if (cmdPipe)
    {
        // Read line number of lines
        for (label cnt = 0; cnt <= line; cnt++)
        {
            char buffer[MAX];
            char* s = fgets(buffer, MAX-1, cmdPipe);
            //cout<< "for line:" << line << " cnt:" << cnt
            //    << " s:" << s << std::endl;
            if (s == NULL)
            {
                return "";
            }
            if (cnt == line)
            {
                string str(buffer);
                return str.substr(0, str.size()-1);
            }            
        }
        pclose(cmdPipe);
    }

    return "";
}


// use popen to call addr2line (using bfd.h directly would have
// meant relinking everything)

void printSourceFileAndLine
(
    Ostream& os,
    const HashTable<label, fileName>& addressMap,
    const fileName& filename,
    const word& address
)
{
    word myAddress = address;

    if (filename.ext() == "so")
    {
        // Convert offset into .so into offset into executable.

        // Check needed since softlinks are resolved in addressMap but
        // not in filename.
        label offset = 0;

        HashTable<label, fileName>::const_iterator fileIter =
            addressMap.find(filename);
        if (fileIter != addressMap.end())
        {
            offset = fileIter();
        }
        else
        {
            HashTable<label, fileName>::const_iterator baseIter =
                addressMap.find(filename.name());
            if (baseIter != addressMap.end())
            {
                offset = baseIter();
            }
        }
        IStringStream addressStr(address.substr(2));
        label addressValue = readHexLabel(addressStr);
        label relativeAddress = addressValue-offset;

        // Reconstruct hex word from address
        OStringStream nStream;
        nStream << "0x" << hex << relativeAddress;
        myAddress = nStream.str();
    }

    if (filename[0]=='/')
    {
        string line = doThePOpenDance
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + myAddress,
            1
        );

        if (line == "??:0")
        {
            //os << " No debug information available";
        }
        else
        {
            string cwdLine(line.replaceAll(cwd()+'/', ""));

            string homeLine(cwdLine.replaceAll(home(), '~'));

            os << " at " << homeLine.c_str();
        }
    }
}


void getSymbolForRaw
(
    Ostream& os,
    const string& raw,
    const fileName& filename,
    const word& address
)
{
    if (filename[0] == '/')
    {
        string fcnt = doThePOpenDance
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + address
        );

        if (fcnt != "")
        {
            os << fcnt.c_str();
            return;
        }
    }
    os << "Uninterpreted: " << raw.c_str();
}

#endif

void error::printStack(Ostream& os)
{
#if defined(HAS_DEMANGLING)

    // Reads the starting addresses for the dynamically linked libraries
    // from the /proc/pid/maps-file
    // I'm afraid this works only for Linux 2.6-Kernels (may work on 2.4)
    // Note2: the filenames in here will have softlinks resolved so will
    // go wrong when having e.g. OpenFOAM installed under a softlink.

    HashTable<label, fileName> addressMap;
    {
        IFstream is("/proc/" + name(pid()) + "/maps");

        while(is.good())
        {
            string line;
            is.getLine(line);

            string::size_type space = line.rfind(' ') + 1;
            fileName libPath = line.substr(space, line.size()-space);

            if (libPath.size() > 0 && libPath[0] == '/')
            {
                string offsetString(line.substr(0,line.find('-')));
                IStringStream offsetStr(offsetString);
                addressMap.insert(libPath, readHexLabel(offsetStr));
            }
        }
    }

    // Get raw stack symbols
    void *array[100];
    size_t size = backtrace(array, 100);
    char **strings = backtrace_symbols(array, size);

    // See if they contain function between () e.g. "(__libc_start_main+0xd0)"
    // and see if cplus_demangle can make sense of part before +
    for (size_t i = 0; i < size; i++)
    {
        string msg(strings[i]);
        fileName programFile;
        word address;

        os << '#' << label(i) << "  ";
        //os << "Raw   : " << msg << "\n\t";
        {
            string::size_type lPos = msg.find('[');
            string::size_type rPos = msg.find(']');
            
            if (lPos != string::npos && rPos != string::npos && lPos<rPos)
            {
                address=msg.substr(lPos+1,rPos-lPos-1);
            }

            string::size_type bracketPos = msg.find('(');
            string::size_type spacePos = msg.find(' ');
            if (bracketPos != string::npos || spacePos != string::npos)
            {
                programFile=msg.substr(0,min(spacePos,bracketPos));

                // not an absolute path
                if (programFile[0] != '/')
                {
                    string tmp = doThePOpenDance("which "+programFile);
                    if (tmp[0]=='/' || tmp[0]=='~')
                    {
                        programFile=tmp;
                    }
                }
            }
        }

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
                    os<< cplusNamePtr;
                    free(cplusNamePtr);
                }
                else
                {
                    os<< cName.c_str();
                }
            }
            else
            {
                string::size_type endBracketPos = msg.find(')', start);
                
                if (endBracketPos != string::size_type(string::npos))
                {
                    string fullName(msg.substr(start, endBracketPos-start));

                    os<< fullName.c_str() << nl;
                }
                else
                {
                    // Print raw message
                    getSymbolForRaw(os, msg, programFile, address);
                }
            }
        }
        else
        {
            // Print raw message
            getSymbolForRaw(os, msg, programFile, address);
        }

        printSourceFileAndLine(os, addressMap, programFile, address);

        os << nl;
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
        printStack(*this);
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
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
        printStack(*this);
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        ::abort();
    }

    if (Pstream::parRun())
    {
        printStack(*this);
        Perr<< endl << *this << endl
            << "\nFOAM parallel run aborting\n" << endl;
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
            printStack(*this);
            Perr<< endl << *this << endl
                << "\nFOAM aborting\n" << endl;
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
