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

#include "clock.H"
#include "string.H"

#include <sstream>
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

char *clock::monthnames[] =
{
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

time_t clock::getTime()
{
    return ::time((time_t*)0);
}


const struct tm clock::rawDate()
{
    time_t t = getTime();
    struct tm *timeStruct = localtime(&t);
    return *timeStruct;
}


string clock::date()
{
    std::ostringstream osBuffer;

    time_t t = getTime();
    struct tm *timeStruct = localtime(&t);

    osBuffer
        << monthnames[timeStruct->tm_mon]
        << ' ' << std::setw(2) << std::setfill('0') << timeStruct->tm_mday
        << ' ' << std::setw(4) << timeStruct->tm_year + 1900;

    return osBuffer.str();
}


string clock::clockTime()
{
    std::ostringstream osBuffer;

    time_t t = getTime();
    struct tm *timeStruct = localtime(&t);

    osBuffer
        << std::setfill('0')
        << std::setw(2) << timeStruct->tm_hour
        << ':' << std::setw(2) << timeStruct->tm_min
        << ':' << std::setw(2) << timeStruct->tm_sec;

    return osBuffer.str();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

clock::clock()
:
    startTime_(getTime()),
    lastTime_(startTime_),
    newTime_(startTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

time_t clock::elapsedClockTime() const
{
    newTime_ = getTime();
    return newTime_ - startTime_;
}


time_t clock::clockTimeIncrement() const
{
    lastTime_ = newTime_;
    newTime_ = getTime();
    return newTime_ - lastTime_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
