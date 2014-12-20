/*----------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2005 OpenCFD Ltd.
    \\/      M anipulation   |
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

Application
    dirToString

Description
    converts a directory path into a string with appropriate capitalisation
    e.g. dir1/dir2 becomes dir1Dir2

Usage
    echo dirName | dirToString

    e.g.
        using csh
        set baseDirName=`echo $dir | sed 's%^\./%%' | $rulesPath/dirToString`

        using ksh
        baseDirName=`echo $dir | sed 's%^\./%%' | $rulesPath/dirToString`

\*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

int main()
{
    int c;
    int nextupper = 0;
    
    while ((c=getchar()) != EOF)
    {
        if (c == '/')
        {
            nextupper = 1;
        }
        else
        {
            if (nextupper)
            {
                putchar(toupper(c));
            }
            else
            {
                putchar(c);
            }
            
            nextupper = 0;
        }
    }

    return 0;
}


/*****************************************************************************/
