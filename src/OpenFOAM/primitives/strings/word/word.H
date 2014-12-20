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
    word

Description
    A class for handling words, derived from string.
    A word is a string of characters containing no white space and
    may be constructed from a string by removing white space.
    Words are delimited be white space.

SourceFiles
    word.C
    wordIO.C

\*---------------------------------------------------------------------------*/

#ifndef word_H
#define word_H

#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class word Declaration
\*---------------------------------------------------------------------------*/

class word
:
    public string
{
    // Private member functions

        //- Strip invalid characters from this word
        void stripInvalid();


public:

    //- Null word
    static const word null;


    // Constructors

        //- Construct null
        inline word();

        //- Construct as copy
        inline word(const word& w);

        //- Construct as copy of character array
        word(const char*);

        //- Construct as copy of string
        word(const string&);

        //- Construct as copy of std::string
        word(const std::string&);

        //- Construct from Istream
        word(Istream& is);


    // Member functions

        //- Is this character valid for a word
        inline static bool valid(char);


    // Member operators

        // Assignment

            void operator=(const word&);
            void operator=(const string&);
            void operator=(const std::string&);
            void operator=(const char*);


    // Friend Operators

        inline friend word operator&(const word&, const word&);


    // IOstream operators

        friend Istream& operator>>(Istream&, word&);
        friend Ostream& operator<<(Ostream&, const word&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
