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

#include "string.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

//- Null string
const string string::null;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Count and return the number of a given character in the string
string::size_type string::count(const char c) const
{
    register size_type cCount=0;

    for
    (
        const_iterator iter = begin();
        iter != end();
        iter++
    )
    {
        if (*iter == c)
        {
            cCount++;
        }
    }

    return cCount;
}


// Replace first occurence of sub-string oldStr with newStr
string& string::replace
(
    const string& oldStr,
    const string& newStr,
    size_type start
)
{
    size_type newStart = start;

    if ((newStart = find(oldStr, newStart)) != npos)
    {
        std::string::replace(newStart, oldStr.size(), newStr);
    }

    return *this;
}


// Replace all occurences of sub-string oldStr with newStr
string& string::replaceAll
(
    const string& oldStr,
    const string& newStr,
    size_type start
)
{
    if (oldStr.size())
    {
        size_type newStart = start;

        while ((newStart = find(oldStr, newStart)) != npos)
        {
            std::string::replace(newStart, oldStr.size(), newStr);
            newStart += newStr.size();
        }
    }

    return *this;
}


// Expand all occurences of environment variables
string& string::expand()
{
    // Expand $VARS

    size_type startEnvar = 0;

    // Repeat until nothing is found
    while
    (
        (startEnvar = find('$', startEnvar)) != npos
     && startEnvar < size()-1
    )
    {
        if (startEnvar == 0 || operator[](startEnvar-1) != '\\')
        {
            // Find end of first occurrence
            size_type endEnvar = npos;
            size_type nd = 0;

            if (operator[](startEnvar+1) == '{')
            {
                endEnvar = find('}', startEnvar);
                nd = 1;
            }
            else
            {
                endEnvar = startEnvar;
                iterator iter = begin() + startEnvar + 1;

                while (iter != end() && (isalnum(*iter) || *iter == '_'))
                {
                    iter++;
                    endEnvar++;
                }
            }

            if (endEnvar != npos && endEnvar != startEnvar)
            {
                string enVar = substr
                (
                    startEnvar + 1 + nd,
                    endEnvar - startEnvar - 2*nd
                );
                string enVarString = getEnv(enVar);

                if (enVarString.size())
                {
                    std::string::replace
                    (
                        startEnvar,
                        endEnvar - startEnvar + 1,
                        enVarString
                    );
                    startEnvar += enVarString.size();
                }
                else
                {
                    startEnvar = endEnvar;
                }
            }
            else
            {
                break;
            }
        }
        else
        {
            startEnvar++;
        }
    }


    // Expand ~ into the home directory path
    string HOME = home();
    startEnvar = 0;

    while ((startEnvar = find('~', startEnvar)) != npos)
    {
        if (startEnvar == 0 || operator[](startEnvar-1) != '\\')
        {
            std::string::replace(startEnvar, 1, HOME);
            startEnvar += HOME.size();
        }
        else
        {
            startEnvar++;
        }
    }

    // Expand initial . into cwd
    if
    (
        size()
     && operator[](0) == '.'
     && (size() == 1 || (size() > 1 && operator[](1) == '/'))
    )
    {
        std::string::replace(0, 1, cwd());
    }

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
