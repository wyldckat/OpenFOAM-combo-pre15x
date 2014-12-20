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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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

    if (controlDict_.found("purgeWrite"))
    {
        purgeWrite_ = readInt(controlDict_.lookup("purgeWrite"));

        if (purgeWrite_ < 0)
        {
            WarningIn("Time::readDict(const dictionary&)")
                << "invalid value for purgeWrite " << purgeWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            purgeWrite_ = 0;
        }

        if (writeControl_ != wcTimeStep && purgeWrite_ > 0)
        {
            FatalIOErrorIn("Time::readDict(const dictionary&)", controlDict_)
                << "writeControl must be set to timeStep for purgeWrite "
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
            WarningIn("Time::readDict(const dictionary&)")
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
        else
        {
            endTime_ = GREAT;
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

        Sout.precision(IOstream::defaultPrecision());
        Serr.precision(IOstream::defaultPrecision());

        Pout.precision(IOstream::defaultPrecision());
        Perr.precision(IOstream::defaultPrecision());
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


bool Time::write
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
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

        timeDict.regIOobject::write
        (
            fmt,
            ver,
            cmp
        );

        bool writeOK = objectRegistry::write
        (
            fmt,
            ver,
            cmp
        );

        if (writeOK && purgeWrite_)
        {
            previousOutputTimes_.push(timeName());

            while(previousOutputTimes_.size() > purgeWrite_)
            {
                rmDir(objectRegistry::path(previousOutputTimes_.pop()));
            }
        }

        return writeOK;
    }
    else
    {
        return false;
    }
}


bool Time::write() const
{
    return regIOobject::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
