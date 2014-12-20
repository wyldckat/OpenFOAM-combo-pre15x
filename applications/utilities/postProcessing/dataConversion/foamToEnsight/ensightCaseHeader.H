if (Pstream::master())
{
    ensightCaseFile << "FORMAT" << nl;
    ensightCaseFile << "type: ensight gold" << nl << nl;

    if (lagrangianExist)
    {
        ensightCaseFile
            << "measured:     1     "
            << prepend + "***." << sprayName
            << nl;
    }
}
