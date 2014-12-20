%{
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
    wmkdep

Description
     A fast dependency list generator which emulates the behaviour and
     output of cpp -M. However, the output contains no duplications and 
     is ~40% faster than cpp.

     The algorithm uses lex to scan for includes and searches the files
     found.  Each file is entererd into a hash table so that files are
     scaned only once.  This is why this program is faster than cpp.

Usage
    wmkdep [ -Idirectory ... -Idirectory] filename

\*----------------------------------------------------------------------------*/

#define FILE_STACK_SIZE 300
#define HASH_TABLE_SIZE 500

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>

void nextFile(const char* fileName);
void importFile(const char* fileName);
void importDirectory(const char* dirName);

#   undef yywrap	/* sometimes a macro by default */

%}

%s NORMAL CMNT CFNAME JFNAME FFNAME
%%

%{
    BEGIN NORMAL;
%}

<NORMAL>"//".*\n      ;                 /* remove c++ style one line comments */

<NORMAL>"/*"          BEGIN CMNT;       /* start removing c style comment */
<CMNT>.               |
<CMNT>\n              ;
<CMNT>"*/"            BEGIN NORMAL;     /* end removing c style comment */

<NORMAL>#[ \t]*include[ \t]+\"  BEGIN CFNAME;
<CFNAME>[^"\n]*       { BEGIN NORMAL; nextFile(yytext); } /*"*/

<NORMAL>^[ \t]*import[ \t]+   BEGIN JFNAME;
<JFNAME>java.*;       { BEGIN NORMAL; }
<JFNAME>org.*;        { BEGIN NORMAL; }
<JFNAME>javax.*;      { BEGIN NORMAL; }
<JFNAME>com.*;        { BEGIN NORMAL; }
<JFNAME>sunw.*;       { BEGIN NORMAL; }
<JFNAME>sun.*;        { BEGIN NORMAL; }
<JFNAME>launcher.*;   { BEGIN NORMAL; }
<JFNAME>[^"\n*]*;     { BEGIN NORMAL; importFile(yytext); } /*"*/
<JFNAME>[^"\n]*\*;    { BEGIN NORMAL; importDirectory(yytext); } /*"*/

<NORMAL>"      "include[ \t]+\'  BEGIN FFNAME;
<FFNAME>[^']*         { BEGIN NORMAL; nextFile(yytext); } /*'*/

<NORMAL>.|\t|\n       ;
%%


int nDirectories;
char** directories;
char* sourceFile = NULL;
char* sourceExt = NULL;
char* objectFile = NULL;
char* classFile = NULL;
char* depFile = NULL;

int main(int argc, char* argv[])
{
    char *dotPtr, *slashPtr;

    if (argc == 1)
    {
        fprintf(stderr,"input file not supplied\n");
        exit(1);
    }

    sourceFile = (char*)malloc(strlen(argv[argc-1]) + 1);
    strcpy(sourceFile, argv[argc-1]);
    fprintf(stderr,"Making dependency list for source file %s\n", sourceFile);
    
    nDirectories = argc-2;
    directories = &argv[1];

    if ((dotPtr = strrchr(sourceFile, '.')) == NULL)
    {
        fprintf
        (
            stderr,
            "Cannot find extension in source file name %s\n",
            sourceFile
        );
        exit(1);
    }

    if ((slashPtr = strrchr(sourceFile, '/')) == NULL)
    {
        slashPtr = sourceFile;
    }
    else
    {
        slashPtr++;
    }

    sourceExt = (char*)malloc(strlen(sourceFile));
    sourceExt[0] = 0;
    strncat
    (
        sourceExt,
        dotPtr+1,
        (&sourceFile[strlen(sourceFile) - 1] - dotPtr)/sizeof(char)
    );

    /* initialise depFile to zero and use strncat rather than strncpy
       because there is a bug in the SGI strncat that if 0 preceeds the .
       it inserts a space */
    depFile = (char*)malloc(strlen(sourceFile) + 3);
    depFile[0] = 0;
    strncat(depFile, sourceFile, (dotPtr - sourceFile)/sizeof(char));
    strcat(depFile, ".dep");

    if (strcmp(sourceExt, "java") == 0)
    {
        classFile = (char*)malloc(strlen(sourceFile) + 17);
        strcpy(classFile, "$(CLASSES_DIR)/");
        strncat(classFile, sourceFile, (dotPtr - sourceFile)/sizeof(char));
        strcat(classFile, ".class");

        printf("%s: %s\n", classFile, depFile);
    }
    else
    {
        objectFile = (char*)malloc(strlen(sourceFile) + 16);
        strcpy(objectFile, "$(OBJECTS_DIR)/");
        strncat(objectFile, slashPtr, (dotPtr - slashPtr)/sizeof(char));
        strcat(objectFile, ".o");

        printf("%s: %s\n", objectFile, depFile);
    }

    nextFile(sourceFile);
    yylex();

    free(sourceFile);
    free(sourceExt);
    free(objectFile);
    free(classFile);
    free(depFile);

    return 0;
}


int currentFile = -1;            /* file pointer stack counter */
FILE* files[FILE_STACK_SIZE];    /* file pointer stack */


/* file name entry in hash table */

struct FileName
{
    char* name;
    struct FileName* next;
};

struct FileName* fileHashTable[HASH_TABLE_SIZE];    /* File hash table */
struct FileName* dirHashTable[HASH_TABLE_SIZE];    /* Directory hash table */


/* lookup name in hash table, if not found insert in table */

int lookUp(struct FileName** hashTable, const char* p)
{
    int ii = 0;
    struct FileName* n;
    struct FileName* nn;

    /* hash */
    const char* pp = p;
    while (*pp) ii = ii<<1 ^ *pp++;
    if (ii < 0) ii = -ii;
    ii %= HASH_TABLE_SIZE;

    /* search */
    for (n=hashTable[ii]; n; n=n->next)
    {
        if (strcmp(p, n->name) == 0)
        {
            /* entry found so return true */
            return 1;
        }
    }

    /* insert */
    nn = (struct FileName*)malloc(sizeof(struct FileName));
    nn->name = (char*)malloc(strlen(p)+1);
    strcpy(nn->name, p);
    nn->next = hashTable[ii];
    hashTable[ii] = nn;

    /* entry not found, and therefore added. return false */
    return 0;
}



/* add a directory name to the file name */

char* addDirectoryName(const char* directoryName, const char* fileName)
{
    char* pathName;
    pathName = (char*)malloc(strlen(directoryName) + strlen(fileName) + 2);

    strcpy(pathName, directoryName);

    if (directoryName[strlen(directoryName)-1] != '/')
    {
        strcat(pathName, "/");
    }

    strcat(pathName, fileName);

    return pathName;
}


/* open a file and put pointer on the file pointer stack */

void nextFile(const char* fileName)
{
    int d;
    char* pathName;

    if (lookUp(fileHashTable, fileName)) return;

    currentFile++;

    if(currentFile == FILE_STACK_SIZE)
    {
        fprintf
        (
            stderr,
            "depth of file search exceeds stack size %d "
            "while opening %s for file %s\n",
            FILE_STACK_SIZE, fileName, sourceFile
        );
        exit(1);
    }

    if(!(files[currentFile] = fopen(fileName, "r")))
    {
        for(d=0; d<nDirectories; d++)
        {
            /* Note the +2 which removes the -I */
            pathName = addDirectoryName(directories[d]+2, fileName);

            if ((files[currentFile] = fopen(pathName, "r")))
            { 
                printf("%s: %s\n", depFile, pathName);
                yyin = files[currentFile];
                free(pathName);

                return;
            }

            free(pathName);
        }
        fprintf
        (
            stderr, 
            "could not open file %s for source file %s\n", 
            fileName, sourceFile
        );

        fflush(stdout);
        currentFile--;
    }
    else
    {
        printf("%s: %s\n", depFile, fileName);
        fflush(stdout);

        yyin = files[currentFile];
    }
}


void dotToSlash(char* fileName)
{
    int i, len;

    len = strlen(fileName);

    for (i=0; i<len; i++)
    {
        if (fileName[i] == '.') fileName[i] = '/';
    }
}


void importFile(const char* fileName)
{
    char* javaFileName;
    int fileNameLen;

    fileNameLen = strlen(fileName);
    javaFileName = (char*)malloc(fileNameLen + 6);
    javaFileName[0] = 0;
    strncat(javaFileName, fileName, fileNameLen - 1);
    dotToSlash(javaFileName);
    strcat(javaFileName, ".java");

    nextFile(javaFileName);

    free(javaFileName);
}


void importDirectory(const char* dirName)
{
    int dirNameLen;
    char *uDirName, *path;
    DIR *source;
    struct dirent *list;

    if (lookUp(dirHashTable, dirName)) return;

    dirNameLen = strlen(dirName);
    uDirName = strdup(dirName);
    dotToSlash(uDirName);
    uDirName[dirNameLen-2] = 0;

    if ((source = opendir(uDirName)) == NULL)
    {
        fprintf
        (
            stderr, 
            "could not open directory %s\n", 
            uDirName
        );
        free(uDirName);

        return;
    }
    else
    {
        /* Read and parse all the entries in the directory */
        while ((list = readdir(source)) != NULL)
        {
            if
            (
                strstr(list->d_name, ".java")
             && !strstr(list->d_name, ".java~")
            )
            {
                path = addDirectoryName(uDirName, list->d_name);
                nextFile(path);
                free(path);
            }
        }

        closedir(source);
        free(uDirName);
    }
}


/* the lexer calls yywrap to handle EOF conditions */

int yywrap()
{
    /* close file which has just reached EOF */
    fclose(files[currentFile]);

    /* set file counter to previous file */
    currentFile--;

    if (currentFile >= 0) /* if file counter refers to a valid file */
    {
        /* reset input file to the previous file on the file pointer stack */
        yyin = files[currentFile];

        /* return 0 to inform lex to continue reading */
        return 0;
    }
    else                 /* else there are no more files on the file stack */
    {
        /* return 1 to inform lex finish now that all files have been read */
        return 1;
    }
}


/*****************************************************************************/
