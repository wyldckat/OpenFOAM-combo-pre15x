#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    Info<< vector::zero << endl
        << vector::one << endl
        << vector::dim << endl
        << vector::rank << endl;

    return 0;
}
