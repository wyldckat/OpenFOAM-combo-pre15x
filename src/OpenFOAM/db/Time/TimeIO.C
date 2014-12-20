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

\*---------------------------------------------------------------------------*/

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Time::readDict()
{
    if (!deltaTchanged_)
    {
        deltaT_ = readScalar(controlDict_.lookup("deltaT"));
    }

    if (controlDict_.found("writeControl"))
    {
        writeControl_ = writeControlNames_.read
        (
            controlDict_.lookup("writeControl")
        );
    }

    if (controlDict_.found("writeInterval"))
    {
        controlDict_.lookup("writeInterval") >> writeInterval_;

        if (writeControl_ == wcTimeStep && label(writeInterval_) < 1)
        {
            FatalIOErrorIn("Time::readDict(const dictionary&)", controlDict_)
                << "writeInterval < 1 for writeControl timeStep"
                << exit(FatalIOError);
        }
    }
    else
    {
        controlDict_.lookup("writeFrequency") >> writeInterval_;
    }

    if (controlDict_.found("cycleWrite"))
    {
        cycleWrite_ = readInt(controlDict_.lookup("cycleWrite"));

        if (cycleWrite_ < 0)
        {
            Warning
                << "Time::readDict(const dictionary&) : "
                << "invalid value for cycleWrite " << cycleWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            cycleWrite_ = 0;
        }

        if (writeControl_ != wcTimeStep && cycleWrite_ > 0)
        {
            FatalIOErrorIn("Time::readDict(const dictionary&)", controlDict_)
                << "writeControl must be set to timeStep for cycleWrite "
                << exit(FatalIOError);
        }
    }

    if (controlDict_.found("timeFormat"))
    {
        word formatName(controlDict_.lookup("timeFormat"));

        if (formatName == "general")
        {
            format_ = general;
        }
        else if (formatName == "fixed")
        {
            format_ = fixed;
        }
        else if (formatName == "scientific")
        {
            format_ = scientific;
        }
        else
        {
            Warning
                << "Time::readDict(const dictionary&) : "
                << "unsupported time format " << formatName
                << endl;
        }
    }

    if (controlDict_.found("timePrecision"))
    {
        precision_ = readLabel(controlDict_.lookup("timePrecision"));
    }

    if (controlDict_.found("stopAt"))
    {
        stopAt_ = stopAtControlNames_.read(controlDict_.lookup("stopAt"));

        if (stopAt_ == saEndTime)
        {
            endTime_ = readScalar(controlDict_.lookup("endTime"));
        }
        else if (stopAt_ == saNoWriteNow)
        {
            endTime_ = value();
        }
        else if (stopAt_ == saWriteNow)
        {
            endTime_ = value();
            outputTime_ = true;
        }
        else if (stopAt_ == saNextWrite)
        {
            endTime_ = GREAT;
        }
        else
        {
            Warning
                << "Time::readDict(const dictionary&) : "
                   "unsupported stopAt option "
                << stopAtControlNames_[stopAt_] << nl
                << "    supported options are "
                   "endTime, noWriteNow, writeNow and nextWrite"
                << nl << endl;
        }
    }
    else
    {
        endTime_ = readScalar(controlDict_.lookup("endTime"));
    }

    dimensionedScalar::name() = timeName(value());

    if (controlDict_.found("writeVersion"))
    {
        writeVersion_ = IOstream::versionNumber
        (
            controlDict_.lookup("writeVersion")
        );
    }

    if (controlDict_.found("writeFormat"))
    {
        writeFormat_ = IOstream::format
        (
            controlDict_.lookup("writeFormat")
        );
    }

    if (controlDict_.found("writePrecision"))
    {
        IOstream::defaultPrecision
        (
            readUint(controlDict_.lookup("writePrecision"))
        );
    }

    if (controlDict_.found("writeCompression"))
    {
        writeCompression_ = IOstream::compression
        (
            controlDict_.lookup("writeCompression")
        );
    }

    if (controlDict_.found("graphFormat"))
    {
        graphFormat_ = word(controlDict_.lookup("graphFormat"));
    }

    if (controlDict_.found("runTimeModifiable"))
    {
        runTimeModifiable_ = Switch(controlDict_.lookup("runTimeModifiable"));
    }
}


bool Time::read()
{
    if (controlDict_.regIOobject::read())
    {
        readDict();
        return true;
    }
    else
    {
        return false;
    }
}


void Time::readModifiedObjects()
{
    if (runTimeModifiable_)
    {
        if (controlDict_.readIfModified())
        {
            readDict();
        }

        objectRegistry::readModifiedObjects();
    }
}


bool Time::write() const
{
    if (outputTime())
    {
        IOdictionary timeDict
        (
            IOobject
            (
                "time",
                timeName(),
                *this
            )
        );

        timeDict.add("index", timeIndex_);
        timeDict.add("deltaT", deltaT_);
        timeDict.add("deltaT0", deltaT0_);

        timeDict.regIOobject::write();

        return regIOobject::write();
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
