#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    Info<< vector::zero << endl
        << vector::one << endl
        << vector::dim << endl
        << vector::rank << endl;

    vector d(0.5, 0.5, 0.5);
    d /= mag(d);

    vector dSmall = (1e-100)*d;
    dSmall /= mag(dSmall);

    Info<< (dSmall - d) << endl;

    return 0;
}
