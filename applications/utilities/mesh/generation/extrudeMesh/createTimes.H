    Info<< nl << "Create databases as time" << endl;

    Time runTime
    (
        Time::controlDictName,
        rootDirSource,
        caseDirSource
    );

    Time runTimeExtruded
    (
        Time::controlDictName,
        rootDirTarget,
        caseDirTarget
    );
