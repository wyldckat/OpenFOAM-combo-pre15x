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

#include <sstream>

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Time, 0);
}


template<>
const char* Foam::NamedEnum<Foam::Time::stopAtControls, 4>::names[] =
{
    "endTime",
    "noWriteNow",
    "writeNow",
    "nextWrite"
};

const Foam::NamedEnum<Foam::Time::stopAtControls, 4>
    Foam::Time::stopAtControlNames_;

template<>
const char* Foam::NamedEnum<Foam::Time::writeControls, 5>::names[] =
{
    "timeStep",
    "runTime",
    "adjustableRunTime",
    "clockTime",
    "cpuTime"
};

const Foam::NamedEnum<Foam::Time::writeControls, 5>
    Foam::Time::writeControlNames_;

Foam::Time::fmtflags Foam::Time::format_(Foam::Time::general);
int Foam::Time::precision_(6);

Foam::word Foam::Time::controlDictName = "controlDict";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Time::adjustDeltaT()
{
    if (writeControl_ == wcAdjustableRunTime)
    {
        scalar timeToNextWrite = max
        (
            0.0,
            (outputTimeIndex_ + 1)*writeInterval_ - (value() - startTime_)
        );

        label nStepsToNextWrite = label(timeToNextWrite/deltaT_ - SMALL) + 1;
        scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

        // Control the increase or decrease of the time step to
        // within a factor of 2.  
        if (newDeltaT >= deltaT_)
        {
            deltaT_ = min(newDeltaT, 2.0*deltaT_);
        }
        else
        {
            deltaT_ = max(newDeltaT, 0.2*deltaT_);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Time::Time
(
    const word& controlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    cycleWrite_(0),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(true)
{
    if (controlDict_.found("startFrom"))
    {
        word startFrom(controlDict_.lookup("startFrom"));

        if (startFrom == "startTime")
        {
            startTime_ = readScalar(controlDict_.lookup("startTime"));
        }
        else
        {
            // Search directory for valid time directories
            instantList Times = findTimes(path());

            if (Times.size() > 0)
            {
                if (startFrom == "firstTime")
                {
                    startTime_ = Times[0].value();
                }
                else if (startFrom == "latestTime")
                {
                    startTime_ = Times[Times.size()-1].value();
                }
                else
                {
                    Warning
                        << "Time::Time(const word&, const fileName&) : " << nl
                        << "    expected startTime, firstTime or latestTime"
                        << " found " << startFrom
                        << " in dictionary " << controlDict_.name() << nl
                        << "    Setting time to 0.0" << endl;
                }
            }
        }
    }
    else if (controlDict_.found("currentTime"))
    {
        token currentTimeToken(controlDict_.lookup("currentTime"));

        if (currentTimeToken.isNumber())
        {
            startTime_ = currentTimeToken.number();
        }
        else if (currentTimeToken.isWord())
        {
            // Search directory for valid time directories
            instantList Times = findTimes(path());

            if (Times.size() > 0)
            {
                if (currentTimeToken.wordToken() == "firstTime")
                {
                    startTime_ = Times[0].value();
                }
                else if (currentTimeToken.wordToken() == "latestTime")
                {
                    startTime_ = Times[Times.size()-1].value();
                }
                else
                {
                    Warning
                        << "Time::Time(const word&, const fileName&) : " << nl
                        << "    expected <time>, firstTime or latestTime"
                        << " found " << currentTimeToken.wordToken()
                        << " in dictionary " << controlDict_.name() << nl
                        << "    Setting time to 0.0" << endl;
                }
            }
        }
    }
    else
    {
        FatalIOError
        (
            "Time::Time(const word&, const fileName&)",
            __FILE__,
            __LINE__,
            controlDict_.name(),
            0
        )   << "    expected startFrom or currentTime (deprecated)"
            << " in dictionary " << controlDict_.name() 
            << exit(FatalIOError);
    }

    setTime(startTime_, 0);

    readDict();
    deltaTSave_ = deltaT_;
    deltaT0_ = deltaTSave_;

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (timeDict.found("deltaT"))
    {
        deltaTSave_ = readScalar(timeDict.lookup("deltaT"));
        deltaT0_ = deltaTSave_;
    }

    if (timeDict.found("index"))
    {
        timeDict.lookup("index") >> timeIndex_;
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::Time::~Time()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::Time::timeName(const scalar t)
{
    std::ostringstream osBuffer;
    osBuffer.setf(ios_base::fmtflags(format_), ios_base::floatfield);
    osBuffer.precision(precision_);
    osBuffer << (1 + SMALL)*t;
    return osBuffer.str();
}


Foam::word Foam::Time::timeName() const
{
    return dimensionedScalar::name();
    //return timeName(timeOutputValue());
}


// Search the construction path for times
Foam::instantList Foam::Time::times() const
{
    return findTimes(path());
}


Foam::word Foam::Time::findInstancePath(const instant& t) const
{
    instantList times = Time::findTimes(rootPath()/caseName());

    forAllReverse(times, i)
    {
        if (times[i] == t)
        {
            return times[i].name();
        }
    }

    return word::null;
}


Foam::instant Foam::Time::findClosestTime(const scalar t) const
{
    instantList times = Time::findTimes(rootPath()/caseName());

    // If there is only one time it's "constant" so return it
    if (times.size() == 1)
    {
        return times[0];
    }

    if (t < times[1].value())
    {
        return times[1];
    }
    else if (t > times[times.size() - 1].value())
    {
        return times[times.size() - 1];
    }

    scalar deltaT = GREAT;
    label closesti = 0;

    for (label i=1; i<times.size(); i++)
    {
        if (mag(times[i].value() - t) < deltaT)
        {
            deltaT = mag(times[i].value() - t);
            closesti = i;
        }
    }

    return times[closesti];
}


Foam::dimensionedScalar Foam::Time::startTime() const
{
    return dimensionedScalar("startTime", dimTime, startTime_);
}


Foam::dimensionedScalar Foam::Time::endTime() const
{
    return dimensionedScalar("endTime", dimTime, endTime_);
}


bool Foam::Time::run() const
{
    return (value() < (endTime_ - 0.5*deltaT_));
}


bool Foam::Time::end() const
{
    return (value() > (endTime_ + 0.5*deltaT_));
}


void Foam::Time::setTime(const Time& t)
{
    value() = t.value();
    dimensionedScalar::name() = t.dimensionedScalar::name();
    timeIndex_ = t.timeIndex_;
}


void Foam::Time::setTime(const instant& inst, const label newIndex)
{
    value() = inst.value();
    dimensionedScalar::name() = inst.name();
    timeIndex_ = newIndex;
}


void Foam::Time::setTime(const dimensionedScalar& newTime, const label newIndex)
{
    setTime(newTime.value(), newIndex);
}


void Foam::Time::setTime(const scalar newTime, const label newIndex)
{
    value() = newTime;

    if (cycleWrite_)
    {
        dimensionedScalar::name() = timeName
        (
            (
                label(startTime_) + timeIndex_/label(writeInterval_) - 1
            )%cycleWrite_ + 1
        );
    }
    else
    {
        dimensionedScalar::name() = timeName(timeToUserTime(newTime));
    }

    timeIndex_ = newIndex;
}


void Foam::Time::setEndTime(const dimensionedScalar& endTime)
{
    setEndTime(endTime.value());
}


void Foam::Time::setEndTime(const scalar endTime)
{
    endTime_ = endTime;
}


void Foam::Time::setDeltaT(const dimensionedScalar& deltaT)
{
    setDeltaT(deltaT.value());
}


void Foam::Time::setDeltaT(const scalar deltaT)
{
    deltaT_ = deltaT;
    deltaTchanged_ = true;
    adjustDeltaT();
}


Foam::TimeState Foam::Time::subCycle(const label nSubCycles)
{
    TimeState ts = *this;

    setTime(*this - deltaT(), timeIndex() - 1);
    deltaT_ /= nSubCycles;
    deltaT0_ /= nSubCycles;

    return ts;
}


void Foam::Time::endSubCycle(const TimeState& ts)
{
    TimeState::operator=(ts);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Time& Foam::Time::operator+=(const dimensionedScalar& deltaT)
{
    readModifiedObjects();

    return operator+=(deltaT.value());
}


Foam::Time& Foam::Time::operator+=(const scalar deltaT)
{
    readModifiedObjects();

    setDeltaT(deltaT);
    operator++();

    return *this;
}


Foam::Time& Foam::Time::operator++()
{
    readModifiedObjects();

    //adjustDeltaT();

    deltaT0_ = deltaTSave_;
    deltaTSave_ = deltaT_;
    setTime(value() + deltaT_, timeIndex_ + 1);

    // If the time is very close to zero reset to zero
    if (mag(value()) < 10*SMALL)
    {
        setTime(0.0, timeIndex_);
    }


    switch(writeControl_)
    {
        case wcTimeStep:
            outputTime_ = !(timeIndex_%label(writeInterval_));
        break;

        case wcRunTime:
        case wcAdjustableRunTime:
        {
            label outputTimeIndex = 
                label(((value() - startTime_) + 0.5*deltaT_)/writeInterval_);

            if (outputTimeIndex > outputTimeIndex_)
            {
                outputTime_ = true;
                outputTimeIndex_ = outputTimeIndex;
            }
            else
            {
                outputTime_ = false;
            }
        }
        break;

        case wcCpuTime:
        {
            label outputTimeIndex =
                label(elapsedCpuTime()/writeInterval_);

            if (outputTimeIndex > outputTimeIndex_)
            {
                outputTime_ = true;
                outputTimeIndex_ = outputTimeIndex;
            }
            else
            {
                outputTime_ = false;
            }
        }
        break;

        case wcClockTime:
        {
            label outputTimeIndex = label(elapsedClockTime()/writeInterval_);
            if (outputTimeIndex > outputTimeIndex_)
            {
                outputTime_ = true;
                outputTimeIndex_ = outputTimeIndex;
            }
            else
            {
                outputTime_ = false;
            }
        }
        break;
    };

    if (stopAt_ == saNextWrite && outputTime_ == true && !end())
    {
        endTime_ = value();
    }

    return *this;
}


Foam::Time& Foam::Time::operator++(int)
{
    return operator++();
}


// ************************************************************************* //
