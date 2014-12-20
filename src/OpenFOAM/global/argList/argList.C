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
    argList

Description
    Make, print and check argList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"
#include "clock.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "foamVersion.H"
#include "JobInfo.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

SLList<string> argList::validArgs;
HashTable<string> argList::validOptions;
HashTable<string> argList::validParOptions;

argList::initValidTables::initValidTables()
{
    validArgs.append("root");
    validArgs.append("case");

    validOptions.insert("parallel", "");
    validParOptions.insert("parallel", "");

    Pstream::addValidParOptions(validParOptions);
}

argList::initValidTables dummyInitValidTables;


int argList::findArg(const string& arg)
{
    int argi = 1;
    for
    (
        SLList<string>::iterator iter = validArgs.begin();
        iter != validArgs.end();
        ++iter, ++argi
    )
    {
        if (iter() == arg)
        {
            return argi;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

argList::argList(int& argc, char**& argv)
:
    args_(argc),
    options_(argc),
    parRun_(false)
{
    wordList execCmpts = fileName(argv[0]).components();
    executable_ = execCmpts[execCmpts.size()-1];

    int nArgs = 0;

    string argListString;

    for (int argi=0; argi<argc; argi++)
    {
        argListString += argv[argi];
        argListString += ' ';

        if (argv[argi][0] == '-')
        {
            const char *optionName = &argv[argi][1];

            if 
            (
                (
                    validOptions.found(optionName)
                 && validOptions[optionName] != ""
                )
             || (
                    validParOptions.found(optionName)
                 && validParOptions[optionName] != ""
                )
            )
            {
                if (argi >= argc - 1)
                {
                    FatalErrorIn("argList::argList(int argc, char* argv[])")
                        << "option " << "'-" << optionName << '\''
                        << " requires an argument"
                        << exit(FatalError);
                }

                options_.insert(optionName, argv[argi + 1]);
                argi++;

                argListString += argv[argi];
                argListString += ' ';
            }
            else
            {
                options_.insert(optionName, "");
            }
        }
        else
        {
            args_[nArgs++] = argv[argi];
        }
    }

    args_.setSize(nArgs);

    Info<< "Exec   : " << argListString.c_str() << endl;


    // For parallel runs prepend Sout output with processor ID
    string procLabel;

    // Case is a single processor run unless it is running parallel
    int nProcs = 1;

    // Check if this run is a parallel run buy searching for any parallel option
    // If found call runPar
    for
    (
        HashTable<string>::iterator iter = options_.begin();
        iter != options_.end();
        ++iter
    )
    {
        if (validParOptions.found(iter.key()))
        {
            parRunControl_.runPar(argc, argv);
            break;
        }
    }

    string dateString = clock::date();
    string timeString = clock::clockTime();
    fileName currentDir = cwd();

    Sout<< procLabel.c_str() << "Date   : " << dateString.c_str() << nl
        << procLabel.c_str() << "Time   : " << timeString.c_str() << nl
        << procLabel.c_str() << "Host   : " << hostName() << nl
        << procLabel.c_str() << "PID    : " << pid() << endl;

    jobInfo.add("startDate", dateString);
    jobInfo.add("startTime", timeString);
    jobInfo.add("userName", userName());
    jobInfo.add("foamVersion", word(FOAMversion));
    jobInfo.add("code", executable_);
    jobInfo.add("argList", argListString);
    jobInfo.add("currentDir", currentDir);
    jobInfo.add("PPID", ppid());
    jobInfo.add("PGID", pgid());

    if (!check())
    {
        FatalError.exit();
    }
    
    // If this actually is a parallel run
    if (parRunControl_.parRun())
    {
        // For the master
        if (Pstream::master())
        {
            fileNameList roots;

            // Set the root and case argument indices
            int rootIndex = findArg("root");
            int caseIndex = findArg("case");

            if (rootIndex && caseIndex)
            {
                rootPath_ = fileName(args_[rootIndex]).expand();;
                globalCase_ = args_[caseIndex];
                case_ = globalCase_;

                IFstream decompostionDictStream
                (
                    rootPath_/globalCase_/"system/decomposeParDict"
                );

                dictionary decompositionDict(decompostionDictStream);

                Switch distributed(false);

                if (decompositionDict.found("distributed"))
                {
                    decompositionDict.lookup("distributed") >> distributed;

                    if (distributed)
                    {
                        decompositionDict.lookup("roots") >> roots;

                        if (roots.size() != Pstream::nProcs())
                        {
                            FatalErrorIn
                                ("argList::argList(int argc, char* argv[])")
                                << "number of entries in "
                                << "decompositionDict::roots"
                                << " is not equal to the number of processors "
                                << Pstream::nProcs()
                                << exit(FatalError);
                        }
                    }
                }

                if (!distributed)
                {
                    label nProcDirs = 0;
                    while 
                    (
                        dir
                        (
                            rootPath_/globalCase_/"processor" 
                          + name(++nProcDirs)
                        )
                    );

                    if (nProcDirs != Pstream::nProcs())
                    {
                        FatalErrorIn("argList::argList(int argc, char* argv[])")
                            << "number of processor directories = "
                            << nProcDirs
                            << " is not equal to the number of processors = "
                            << Pstream::nProcs()
                            << exit(FatalError);
                    }
                }
            }

            // Distribute the master's argument list
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                if (rootIndex && caseIndex && roots.size())
                {
                    args_[rootIndex] =  roots[slave-1];
                }

                OPstream toSlave(slave, 0, false);
                toSlave << args_ << options_;
            }

            if (rootIndex && caseIndex && roots.size())
            {
                args_[rootIndex] = rootPath_;
            }
        }
        else
        {
            // Collect the master's argument list
            IPstream fromMaster(Pstream::masterNo());
            fromMaster >> args_ >> options_;
        }

        procLabel = '[' + word(name(Pstream::myProcNo())) + "] ";
        nProcs = Pstream::nProcs();
    }

    // Do transfer of machineInfo from slave to master.
    wordList machinesInfo(1);
    machinesInfo[0] = hostName();

    labelList pids(1);
    pids[0] = pid();

    if (parRunControl_.parRun())
    {
        if (Pstream::master())
        {
            machinesInfo.setSize(Pstream::nProcs());
            pids.setSize(Pstream::nProcs());

            // From all slaves
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                label slavei = slave - Pstream::masterNo();
                IPstream fromSlave(slave);
                fromSlave >> machinesInfo[slavei] >> pids[slavei];
            }
        }
        else
        {
            OPstream toMaster(Pstream::masterNo(), 0, false);
            toMaster << machinesInfo[0] << pids[0];
        }
    }

    wordList slaveProcs(machinesInfo.size() - 1);
    for (int proci = 1; proci < machinesInfo.size(); proci++)
    {
         slaveProcs[proci-1] = machinesInfo[proci] + "." + name(pids[proci]);
    }

    if (int rootIndex = findArg("root"))
    {
        rootPath_ = fileName(args_[rootIndex]).expand();
    }

    if (int caseIndex = findArg("case"))
    {
        globalCase_ = args_[caseIndex];

        if (globalCase_[globalCase_.size() - 1] == '/')
        {
            FatalErrorIn("argList::argList(int& argc, char**& argv)")
                << "case " << globalCase_
                << " does not exist in root " << rootPath_
                << exit(FatalError);
        }

        if (parRunControl_.parRun())
        {
            case_ =
                globalCase_/::Foam::word(::Foam::word("processor")
              + name(Pstream::myProcNo()));
        }
        else
        {
            case_ = globalCase_;
        }
    }


    Sout<< procLabel.c_str() << "Root   : " << rootPath_.c_str() << nl
        << procLabel.c_str() << "Case   : " << globalCase_.c_str() << nl
        << procLabel.c_str() << "Nprocs : " << nProcs << endl;

    if (slaveProcs.size() != 0)
    {
        Sout<< procLabel.c_str() << "Slaves : " << slaveProcs << endl;
    }

    jobInfo.add("root", rootPath_);
    jobInfo.add("case", globalCase_);
    jobInfo.add("nProcs", nProcs);
    if (slaveProcs.size())
    {
        jobInfo.add("slaves", slaveProcs);
    }
    jobInfo.write();

    // Set the root and case as environment variables
    setEnv("FOAM_ROOT", rootPath_, true);
    setEnv("FOAM_CASE", globalCase_, true);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

argList::~argList()
{
    jobInfo.end();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void argList::noParallel()
{
    validOptions.erase("parallel");
}


void argList::printUsage() const
{
    Info<< "Usage: " << executable_;

    for
    (
        SLList<string>::iterator iter = validArgs.begin();
        iter != validArgs.end();
        ++iter
    )
    {
        Info<< ' ' << '<' << iter().c_str() << '>';
    }

    for
    (
        HashTable<string>::iterator iter = validOptions.begin();
        iter != validOptions.end();
        ++iter
    )
    {
        Info<< " [-" << iter.key();

        if (iter().size())
        {
            Info<< ' ' << iter().c_str();
        }

        Info<< ']';
    }

    Info<< endl;
}


bool argList::check() const
{
    if (Pstream::master())
    {
        bool ok = true;

        if (args_.size() - 1 != validArgs.size())
        {
            FatalErrorIn("argList::check() const")
                << "Wrong number of arguments, expected " << validArgs.size()
                << " found " << args_.size() - 1 << endl;
            ok = false;
        }
        else
        {
            forAllConstIter(HashTable<string>, options_, iter)
            {
                if
                (
                    !validOptions.found(iter.key())
                 && !validParOptions.found(iter.key())
                )
                {
                    FatalErrorIn("argList::check() const")
                        << "Invalid option: -" << iter.key() << endl;
                    ok = false;
                }
            }
        }

        if (!ok)
        {
            Info<< endl;
            printUsage();
        }

        return ok;
    }
    else
    {
        return true;
    }
}


bool argList::checkRootCase() const
{
    if (args_.size() >= 2)
    {
        if (!dir(rootPath()))
        {
            FatalErrorIn("argList::checkRootCase() const")
                << executable_
                << ": cannot open root directory " << rootPath()
                << endl;

            return false;
        }

        if (!dir(path()))
        {
            FatalErrorIn("argList::checkRootCase() const")
                << executable_
                << ": Cannot open case directory " << path()
                << endl;

            return false;
        }

        return true;
    }
    else
    {
        FatalErrorIn("argList::checkRootCase() const")
            << executable_
            << ": <root> <case> not given" << endl;

        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
