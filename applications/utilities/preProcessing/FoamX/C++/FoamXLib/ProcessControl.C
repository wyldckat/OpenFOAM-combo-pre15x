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
// Standard header files.
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <errno.h>

// Foam header files.
#include "ProcessControl.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "timer.H"
#include "long.H"

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXNameSpaces.H"


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void FoamX::ProcessControl::readEntries(const Foam::dictionary& procDict)
{
    if (procDict.found("remoteShell"))
    {
        remoteShell_ = Foam::string(procDict.lookup("remoteShell"));
    }

    if (procDict.found("remoteCp"))
    {
        remoteShell_ = Foam::string(procDict.lookup("remoteCp"));
    }

    if (procDict.found("timeOut"))
    {
        timeOut_ = readLabel(procDict.lookup("timeOut"));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
FoamX::ProcessControl::ProcessControl
(
    const Foam::fileName& fxSystemConfigFileName,
    const Foam::fileName& fxUserConfigFileName
)
:
    remoteShell_("rsh"),
    remoteCp_("rcp"),
    timeOut_(60)
{
    static const char* functionName = 
        "FoamX::ProcessControl::ProcessControl"
        "(const Foam::fileName& fxSystemConfigFileName)";

    if (!Foam::exists(fxSystemConfigFileName))
    {
        throw FoamXError
        (
            E_FAIL,
            "Configuration dictionary '" + fxSystemConfigFileName
          + "' not found",
            functionName,
            __FILE__, __LINE__
        );
    }

    Foam::dictionary configDict((Foam::IFstream(fxSystemConfigFileName)()));
    readEntries(configDict.subDict("processControl"));

    if (Foam::exists(fxUserConfigFileName))
    {
        Foam::dictionary configDict((Foam::IFstream(fxUserConfigFileName)()));

        if (configDict.found("processControl"))
        {
            readEntries(configDict.subDict("processControl"));
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::stringList FoamX::ProcessControl::remoteShellArgs
(
    const Foam::string& userName,
    const Foam::string& hostName,
    const Foam::stringList& arguments,
    const Foam::string& logName,
    const bool backGround
) const
{
    // Get length of resulting argument string

    unsigned int argsLen = arguments.size();
    if (backGround)
    {
        argsLen++;
    }

    bool remote = (hostName != Foam::hostName());
    if (remote)
    {
        argsLen += 4;
    }

    bool haveLogName = false;
    if (logName.length() != 0)
    {
        haveLogName = true;

        argsLen+= 2;
    }

    // Now we have collected
    //  - argsLen
    //  - remote        (true if on remote machine)
    //  - haveLogName   (true if output needs to go into log)


    Foam::stringList args(argsLen);

    // Remote invocation
    unsigned int argi = 0;
    if (remote)
    {
        args[argi++] = remoteShell_;
        args[argi++] = hostName;
        args[argi++] = "-l";
        args[argi++] = userName;
    }

    // Arguments
    for (int i = 0; i <arguments.size(); i++)
    {
        args[argi++] = arguments[i];
    }

    // Redirection
    if (haveLogName)
    {
        // Don't redirect stderr.
        if (remote)
        {
            // Quoted ">" so redirection is done remotely
            args[argi++] = "\">\"";
            args[argi++] = logName;
        }
        else
        {
            args[argi++] = ">";
            args[argi++] = logName;
        }
    }

    if (backGround)
    {
        args[argi++] = "&";
    }
    return args;
}


Foam::stringList FoamX::ProcessControl::remoteCpArgs
(
    const Foam::string& userName,
    const Foam::string& hostName,
    const Foam::fileName& src,
    const Foam::fileName& dest
) const
{
    if (hostName == Foam::hostName())
    {
        Foam::stringList args(4);
        args[0] = "cp";
        args[1] = "-r";
        args[2] = src;
        args[3] = dest;

        return args;
    }
    else
    {
        Foam::stringList args(4);
        args[0] = remoteCp_;
        args[1] = "-r";
        args[2] = userName + '@' + hostName + ':' + src;
        args[3] = dest;

        return args;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
       
// Return the command line as a string
Foam::string FoamX::ProcessControl::commandString
(
    const Foam::stringList& argList
)
{
    Foam::string cms = argList[0];

    for (label i=1; i<argList.size(); i++)
    {
        cms += ' ' + argList[i];
    }

    return cms;
}


// Wait for process to finish.
pid_t FoamX::ProcessControl::waitpid(pid_t pid)
{
    return ::waitpid(pid, NULL, 0);
}


pid_t FoamX::ProcessControl::fork
(
    const Foam::stringList& argList,
    const Foam::fileName& logFile
)
{
    pid_t pid = -1;

    // We are in the forked process so execute the command.
    // Allocate a buffer to hold the command arguments.
    List<char*> argv(argList.size() + 1);

    forAll(argList, i)
    {
        argv[i] = (char*)argList[i].c_str();
    }
    argv[argList.size()] = NULL;

    // Fork this process.
    if ((pid = ::fork()) < 0)
    {
        return -1;
    }

    // Check that the fork worked.
    if (pid == 0)
    {
        if (FoamXError::debug)
        {
            if (logFile.size() != 0)
            {
                Info<< "Running Process [" << getpid() << "] "
                << commandString(argList) << " > " << logFile << endl;
            }
            else
            {
                Info<< "Running Process [" << getpid() << "] "
                    << commandString(argList) << endl;
            }
        }

        if (logFile.size() != 0)
        {
            int outFd =
                open
                (
                    logFile.c_str(),
                    O_WRONLY | O_CREAT | O_TRUNC,
                    S_IRWXU | S_IRWXG | S_IRWXO
                );

            if (dup2(outFd, STDOUT_FILENO) != STDOUT_FILENO)
            {
                perror("ProcessControl::fork : dup2 error to stdout");
            }
            close (outFd);
        }

        // Replace this process with the required command.
        execvp(argv[0], argv.begin());

        // execvp failed so exit this process and return a -1 error code.
        ::_exit(-1);
    }

    return pid;
}


int FoamX::ProcessControl::system
(
    const Foam::stringList& argList,
    const int timeOut
)
{
    static const char* functionName = 
        "FoamX::ProcessControl::system"
        "(const Foam::stringList&, const int)";

    Foam::string cmd(commandString(argList));

    Info<< "Doing (with timeout " << timeOut << ") : " << cmd << endl;

    if (timeOut <= 0)
    {
        return ::system(cmd.c_str());
    }
    else
    {
        char *argp[] =
        {
            "/bin/sh",
            "-c",
            (char *)cmd.c_str(),
            NULL
        };

        pid_t pid = -1;

        if ((pid = ::fork()) < 0)
        {
            if (FoamXError::debug)
            {
                Info<< functionName
                    << " fork failed." << endl;
            }
            return -1;
        }

        if (pid == 0)
        {
            execv("/bin/sh", argp);

            // Should never come here.
            _exit(127);
            return -1;
        }
        else
        {
            if (FoamXError::debug)
            {
                Info<< functionName
                    << " forked pid:" << pid << "." << endl;
            }

            timer myTimer(timeOut);

            if (timedOut(myTimer))
            {
                // Timed out
                FoamX::ProcessControl::kill(pid, SIGTERM);
                return TIMEDOUT;
            }

            int pStat;
            while (::waitpid(pid, &pStat, 0) == -1)
            {
                if (errno != EINTR) 
                {
                    pStat = -1;
                    break;
                }
            }

            //if (FoamXError::debug)
            //{
                Info<< "Finished doing (with timeout " << timeOut << ") : "
                    << cmd << endl;
            //}

            return pStat;
        }
    }
}



int FoamX::ProcessControl::kill(pid_t pid, int sig)
{
    int pStat = ::kill(pid, sig);

    timer myTimer(60);

    if (timedOut(myTimer))
    {
        return -1;
    }
    
    if (::waitpid(pid, &pStat, 0) == -1)
    {
        return pStat;
    }
    else
    {
        return 0;
    }
}


int FoamX::ProcessControl::kill(const word& host, pid_t pid, int sig) const
{
    if (host == Foam::hostName())
    {
        return FoamX::ProcessControl::kill(pid, sig);
    }
    else
    {
        stringList args(3);
        args[0] = "kill";
        args[1] = "-" + Foam::name(sig);
        args[2] = Foam::name(pid);

        stringList remoteArgs
        (
            remoteShellArgs
            (
                Foam::userName(),
                host,
                args,
                "",
                false
            )
        );

        return system
        (
            remoteArgs,
            timeOut()
        );
    }
}


int FoamX::ProcessControl::suspend(const word& host, pid_t pid) const
{
    return kill(host, pid, SIGSTOP);
}


int FoamX::ProcessControl::cont(const word& host, pid_t pid) const
{
    return kill(host, pid, SIGCONT);
}


// ************************************************************************* //
