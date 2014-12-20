HashTable<bool> volScalarHash;
HashTable<bool> volVectorHash;
HashTable<bool> surfScalarHash;
HashTable<bool> surfVectorHash;
HashTable<bool> sprayScalarHash;
HashTable<bool> sprayVectorHash;

forAll(Times, timeI)
{
    runTime.setTime(Times[timeI], timeI);

    // Add all fields to hashtable
    IOobjectList objects(mesh, runTime.timeName());
    {
        wordList fieldNames(objects.names(volScalarField::typeName));
        forAll(fieldNames, fieldI)
        {
            if (!volScalarHash.found(fieldNames[fieldI]))
            {
                volScalarHash.insert(fieldNames[fieldI], true);
            }
        }
    }
    {
        wordList fieldNames(objects.names(volVectorField::typeName));
        forAll(fieldNames, fieldI)
        {
            if (!volVectorHash.found(fieldNames[fieldI]))
            {
                volVectorHash.insert(fieldNames[fieldI], true);
            }
        }
    }
    {
        wordList fieldNames(objects.names(surfaceScalarField::typeName));
        forAll(fieldNames, fieldI)
        {
            if (!surfScalarHash.found(fieldNames[fieldI]))
            {
                surfScalarHash.insert(fieldNames[fieldI], true);
            }
        }
    }
    {
        wordList fieldNames(objects.names(surfaceVectorField::typeName));
        forAll(fieldNames, fieldI)
        {
            if (!surfVectorHash.found(fieldNames[fieldI]))
            {
                surfVectorHash.insert(fieldNames[fieldI], true);
            }
        }
    }


    // Same for spray
    IOobjectList sprayObjects(mesh, runTime.timeName(), "lagrangian");
    {
        wordList fieldNames(sprayObjects.names(scalarIOField::typeName));
        forAll(fieldNames, fieldI)
        {
            if (!sprayScalarHash.found(fieldNames[fieldI]))
            {
                sprayScalarHash.insert(fieldNames[fieldI], true);
            }
        }
    }
    {
        wordList fieldNames(sprayObjects.names(vectorIOField::typeName));
        forAll(fieldNames, fieldI)
        {
            if (!sprayVectorHash.found(fieldNames[fieldI]))
            {
                sprayVectorHash.insert(fieldNames[fieldI], true);
            }
        }
    }
}


wordList volScalarNames(volScalarHash.toc());
wordList volVectorNames(volVectorHash.toc());
wordList surfScalarNames(surfScalarHash.toc());
wordList surfVectorNames(surfVectorHash.toc());
wordList sprayScalarNames(sprayScalarHash.toc());
wordList sprayVectorNames(sprayVectorHash.toc());
