    scalar newTotalVolume = sum(mesh.V());
    scalar totalVolRatio = newTotalVolume/totalVolume;

    Info << "Total volume error: " << totalVolRatio - 1 << endl;

    totalVolume = newTotalVolume;
