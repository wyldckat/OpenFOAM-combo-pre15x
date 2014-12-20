#include "primitiveFields.H"
#include "cpuTime.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    const label size = 10000;

    double f1[size], f2[size], f3[size], f4[size];

    for (register int i=0; i<size; i++)
    {
        f1[i] = 1.0;
        f2[i] = 1.0;
        f3[i] = 1.0;
    }

    cpuTime executionTime1;

    for (int j=0; j<10000; j++)
    {
        for (register int i=0; i<size; i++)
        {
            f4[i] = f1[i] + f2[i] - f3[i];
        }
    }

    Info<< "ExecutionTime = "
        << executionTime1.elapsedCpuTime()
        << " s\n" << endl << endl;

    Info << f4[1] << endl;


    scalarField sf1(size, 1.0), sf2(size, 1.0), sf3(size, 1.0), sf4(size);

    cpuTime executionTime2;

    for (register int j=0; j<10000; j++)
    {
        sf4 = sf1 + sf2 - sf3;
    }

    Info<< "ExecutionTime = "
        << executionTime2.elapsedCpuTime()
        << " s\n" << endl << endl;

    Info << sf4[1] << endl;
}
