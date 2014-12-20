#include "fvCFD.H"

word itoa(const label n)
{
    const label offset = '0';
    const label length = 3;
    char val[length];
    label leftOfN = n;
    for(label i=0;i<length;i++)
    {
        label j = (label)leftOfN/pow(10,length-i-1);
        leftOfN -= j*pow(10,length-i-1);
        val[i] = offset + j;
    }
    val[length] = 0;
    return val;
}
