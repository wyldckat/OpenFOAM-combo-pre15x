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
    ISstream

Description
    Generic input stream.

SourceFiles
    ISstreamI.H
    ISread.C
    ISreadToken.C

\*---------------------------------------------------------------------------*/

#ifndef ISstream_H
#define ISstream_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ISstream Declaration
\*---------------------------------------------------------------------------*/

class ISstream
:
    public Istream
{
    // Private data

        string name_;
        istream& is_;


    // Private member functions

        char nextValid();


protected:

    // Protected member functions

        //- Return the stream reference
        istream& stream()
        {
            return is_;
        }

        //- Return the const stream reference
        const istream& stream() const
        {
            return is_;
        }


public:

    // Constructors

        //- Construct as wrapper around istream
        ISstream
        (
            istream& is,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        )
        :
            Istream(format, version, compression),
            name_(name),
            is_(is)
        {
            if (is_.good())
            {
                setOpened();
                setGood();
            }
            else
            {
                setState(is_.rdstate());
            }
        }


    // Destructor

        virtual ~ISstream()
        {}


    // Member functions

        // Inquirery

            //- Return the name of the stream
            //  Useful for Fstream to return the filename
            virtual const string& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            //  Useful to alter the stream name
            virtual string& name()
            {
                return name_;
            }

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const;


        // Read functions

            //- Raw, low-level get character function.
            inline ISstream& get(char&);

            //- Raw, low-level putback character function.
            inline ISstream& putback(const char&);

            //- Return next token from stream
            virtual Istream& read(token&);

            //- Read a character
            virtual Istream& read(char&);

            //- Read a word
            virtual Istream& read(word&);

            // Read a string (including enclosing double-quotes)
            virtual Istream& read(string&);

            //- Read a label
            virtual Istream& read(label&);

            //- Read a scalar
            virtual Istream& read(scalar&);

            //- Read binary block
            virtual Istream& read(char*, std::streamsize);

            //- Rewind and return the stream so that it may be read again
            virtual Istream& rewind();


        // Stream state functions

            //- Set flags of output stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags flags);


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
