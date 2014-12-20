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
    triSurface

Description
    
\*---------------------------------------------------------------------------*/

%{

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include "IStringStream.H"
#include "triSurface.H"
#include "STLpoint.H"
#include "SLList.H"
#include "fileName.H"
#include "hashSignedLabel.H"
#include "wordList.H"

using namespace Foam;

label nTriangles = 0;
short region = -1;      /* current region */
short maxRegion = 0;    /* max region */

label lineNo = 1;
word startError;

SLList<STLpoint> STLpoints;
SLList<STLpoint> STLnormals;
SLList<label> STLlabels;
HashTable<label, word> STLsolidNames;

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
extern "C" int yywrap()
{
    return 1;
}

%}

one_space             [ \t\f\015]
space                 {one_space}*
some_space            {one_space}+
cspace                ","{space}

alpha                 [_A-Za-z]
digit                 [0-9]
dec_digit             [0-9]
octal_digit           [0-7]
hex_digit             [0-9a-fA-F]

identifier            {alpha}({alpha}|{digit})*
integer               {dec_digit}+
label                 [1-9]{dec_digit}*
zeroLabel             {digit}*

word                  ([[:alnum:]]|[[:punct:]])*
string                {word}({some_space}{word})*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

double                (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))
float                 {double}

x                     {float}
y                     {float}
z                     {float}

solid                 {space}("solid"|"SOLID"){space}
color                 {space}("color"|"COLOR"){some_space}{float}{some_space}{float}{some_space}{float}{space}
facet                 {space}("facet"|"FACET"){space}
normal                {space}("normal"|"NORMAL"){space}
point                 {space}{x}{some_space}{y}{some_space}{z}{space}
outerloop             {space}("outer"{some_space}"loop")|("OUTER"{some_space}"LOOP"){space}
vertex                {space}("vertex"|"VERTEX"){space}
endloop               {space}("endloop"|"ENDLOOP"){space}
endfacet              {space}("endfacet"|"ENDFACET"){space}
endsolid              {space}("endsolid"|"ENDSOLID")({some_space}{word})*


 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%option stack

%x readSolidName
%x readFacet
%x readNormal
%x readVertices
%x readVertex
%x stlerror

%%

%{
   
STLpoint normal;

static const char* stateNames[7] =
{
    "reading solid",
    "reading solid name",
    "reading facet",
    "reading normal",
    "reading vertices",
    "reading vertex",
    "error"
};

static const char* stateExpects[7] =
{
    "'solid', 'color', 'facet' or 'endsolid'",
    "<string>",
    "'normal', 'outer loop' or 'endfacet'",
    "<float> <float> <float>",
    "'vertex' or 'endloop'",
    "<float> <float> <float>",
    ""
};


%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

 /*                      ------ Reading control header ------                 */

{solid} {
        BEGIN(readSolidName);
    }

<readSolidName>{string} {
        word solidName(YYText());
        if (STLsolidNames.found(solidName))
        {
            region = STLsolidNames[solidName];
        }
        else
        {
            region = maxRegion++;
            STLsolidNames.insert(solidName, region);
        }
        BEGIN(INITIAL);
    }

{color} {
    }

{facet} {
        BEGIN(readFacet);
    }

<readFacet>{normal} {
        BEGIN(readNormal);
    }

<readNormal>{point} {
        IStringStream normalStream(YYText());
        normalStream >> normal;
        STLnormals.append(normal);
        BEGIN(readFacet);
    }

<readFacet>{outerloop} {
        BEGIN(readVertices);
    }

<readVertices>{vertex} {
        BEGIN(readVertex);
    }

<readVertex>{point} {
        IStringStream pointStream(YYText());
        STLpoints.append(STLpoint(pointStream));
        BEGIN(readVertices);
    }

<readVertices>{endloop} {
        BEGIN(readFacet);
    }

<readFacet>{endfacet} {
        nTriangles++;
        STLlabels.append(region);
        BEGIN(INITIAL);
    }

{endsolid} {
    }


 /* ------------------ Ignore remaining space and \n s. --------------------- */

<*>{space} {}
<*>\n      { lineNo++; }

 /* ------ Ignore remaining space and \n s.  Any other characters are errors. */

<*>. {
        startError = YYText();
        yy_push_state(stlerror);
    }

<stlerror>.* {
        yy_pop_state();
        FatalErrorIn
        (
            "triSurface::readSTLASCII(const fileName& STLfileName)"
        )   << "while " << stateNames[YY_START] << " on line " << lineNo << nl
            << "    expected " << stateExpects[YY_START]
            << " but found '" << startError.c_str() << YYText() << "'"
            << exit(FatalError);
    }


 /*  ------------------------ On EOF terminate. ----------------------------  */

<<EOF>> {
            yyterminate();
    }
%%


#include <fstream>

bool triSurface::readSTLASCII(const fileName& STLfileName)
{
    std::ifstream STLstream(STLfileName.c_str());

    if (!STLstream)
    {
        FatalErrorIn
        (
            "triSurface::readSTLASCII(const fileName& STLfileName)"
        )   << "file " << STLfileName << " not found"
            << exit(FatalError);
    }

    nTriangles = 0;
    region = -1;      /* current region */
    maxRegion = 0;    /* max region */

    lineNo = 1;

    STLpoints.clear();
    STLnormals.clear();
    STLlabels.clear();
    STLsolidNames.clear();


    yyFlexLexer lexer(&STLstream);
    while(lexer.yylex() != 0);

    if (STLpoints.size() != 3*STLnormals.size())
    {
        FatalErrorIn
        (
            "triSurface::readSTLASCII(const fileName& STLfileName)"
        )   << "in file " << STLfileName << endl
            << "Problem: read " << STLnormals.size() << " normals"
            << " but " << STLpoints.size() << " points"
            << exit(FatalError);
    }

    pointField rawPoints(STLpoints.size());

    label i = 0;
    for
    (
        SLList<STLpoint>::iterator iter = STLpoints.begin();
        iter != STLpoints.end();
        ++iter
    )
    {
        rawPoints[i++] = *iter;
    }

    pointField rawNormals(STLnormals.size());
    i = 0;
    for
    (
        SLList<STLpoint>::iterator iter = STLnormals.begin();
        iter != STLnormals.end();
        ++iter
    )
    {
        rawNormals[i++] = *iter;
    }

    setSize(nTriangles);

    label pointI = 0;
    SLList<label>::iterator iter = STLlabels.begin();
    forAll (*this, i)
    {
        operator[](i)[0] = pointI++;
        operator[](i)[1] = pointI++;
        operator[](i)[2] = pointI++;
        operator[](i).region() = *iter;
        ++iter;
    }

    stitchTriangles(rawPoints);

    // Convert solidNames into regionNames
    patches_.setSize(STLsolidNames.size());

    for
    (
        HashTable<label, word>::const_iterator iter = STLsolidNames.begin();
        iter != STLsolidNames.end();
        ++iter
    )
    {
        patches_[iter()].name() = iter.key();
    }

    setDefaultPatches();

    return true;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of STLToFoam.L
 \* ------------------------------------------------------------------------- */
