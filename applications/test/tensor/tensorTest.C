#include "tensor.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    tensor t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    tensor t2(1, 2, 3, 1, 2, 3, 1, 2, 3);

    tensor t3 = t1 + t2;

    Info<< t3 << endl;

    tensor t4(3,-2,1,-2,2,0,1, 0, 4);

    Info<< hinv(t4) << endl;
    Info<< (hinv(t4) & t4) << endl;

    Info<< t1.x() << t1.y() << t1.z() << endl;

    return(0);
}
