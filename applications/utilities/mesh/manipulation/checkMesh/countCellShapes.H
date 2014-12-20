    // Construct shape recognizers
    hexMatcher hex;
    prismMatcher prism;
    wedgeMatcher wedge;
    pyrMatcher pyr;
    tetWedgeMatcher tetWedge;
    tetMatcher tet;

    // Counters for different cell types
    label nHex = 0;
    label nWedge = 0;
    label nPrism = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTetWedge = 0;
    label nUnknown = 0;

    for(label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (hex.isA(mesh, cellI))
        {
            nHex++;
        }
        else if (tet.isA(mesh, cellI))
        {
            nTet++;
        }
        else if (pyr.isA(mesh, cellI))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, cellI))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, cellI))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, cellI))
        {
            nTetWedge++;
        }
        else
        {   
            nUnknown++;
        }
    }

    Pout<< "Number of cells by type: " << nl
        << "    hexahedra:     " << nHex << nl
        << "    prisms:        " << nPrism << nl
        << "    wedges:        " << nWedge << nl
        << "    pyramids:      " << nPyr << nl
        << "    tet wedges:    " << nTetWedge << nl
        << "    tetrahedra:    " << nTet << nl
        << "    polyhedra:     " << nUnknown
        << endl;

