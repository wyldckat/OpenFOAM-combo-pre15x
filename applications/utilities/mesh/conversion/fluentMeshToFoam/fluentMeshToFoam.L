/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    fluentMeshToFoam

Description
	Converts a Fluent mesh to FOAM format. 

\*---------------------------------------------------------------------------*/

%{

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include "argList.H"
#include "Time.H"
#include "IStringStream.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellShape.H"

#include "readHexLabel.H"
#include "cellShapeRecognition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

const scalar convertToMeters = 1.0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label dimensionOfGrid = 0;

label nPoints = 0;
label nFaces = 0;
label nCells = 0;

pointField points(0);
SLList<label> pointGroupZoneID;
SLList<label> pointGroupStartIndex;
SLList<label> pointGroupEndIndex;


faceList faces(0);
labelList owner(0);
labelList neighbour(0);
SLList<label> faceGroupZoneID;
SLList<label> faceGroupStartIndex;
SLList<label> faceGroupEndIndex;

labelList fluentCellModelID(0);
SLList<label> cellGroupZoneID;
SLList<label> cellGroupStartIndex;
SLList<label> cellGroupEndIndex;
SLList<label> cellGroupType;

// number of zones adjusted at run-time if necessary
label maxZoneID = 100; 
label zoneIDBuffer = 10;

wordList patchTypeIDs(maxZoneID);
wordList patchNameIDs(maxZoneID);

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
extern "C" int yywrap()
{
    return 1;
}

%}

one_space                  [ \t\f\015]
space                      {one_space}*
some_space                 {one_space}+
cspace                     ","{space}
spaceNl                    ({space}|\n)*

alpha                      [_[:alpha:]]
digit                      [[:digit:]]
decDigit                   [[:digit:]]
octalDigit                 [0-7]
hexDigit                   [[:xdigit:]]

lbrac                      "("
rbrac                      ")"
quote                      \"
dash                       "-"
dotColonDash               [.:-]

schemeSpecialInitial       [!$%&*/:<=>?~_^#.]
schemeSpecialSubsequent    [.+-]
schemeSymbol               (({some_space}|{alpha}|{quote}|{schemeSpecialInitial})({alpha}|{quote}|{digit}|{schemeSpecialInitial}|{schemeSpecialSubsequent})*)


identifier                 {alpha}({alpha}|{digit})*
integer                    {decDigit}+
label                      [1-9]{decDigit}*
hexLabel                   {hexDigit}+
zeroLabel                  {digit}*

word                       ({alpha}|{digit}|{dotColonDash})*

exponent_part              [eE][-+]?{digit}+
fractional_constant        [-+]?(({digit}*"."{digit}+)|({digit}+".")|({digit}))

double                     ((({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))|0)

x                          {double}
y                          {double}
z                          {double}
scalar                     {double}
labelListElement           {space}{zeroLabel}
hexLabelListElement        {space}{hexLabel}
scalarListElement          {space}{double}
schemeSymbolListElement    {space}{schemeSymbol}
labelList                  ({labelListElement}+{space})
hexLabelList               ({hexLabelListElement}+{space})
scalarList                 ({scalarListElement}+{space})
schemeSymbolList           ({schemeSymbolListElement}+{space})

starStar                   ("**")
text                       ({space}({word}*{space})*)

dateDDMMYYYY               ({digit}{digit}"/"{digit}{digit}"/"{digit}{digit}{digit}{digit})
dateDDMonYYYY              ((({digit}{digit}{space})|({digit}{space})){alpha}*{space}{digit}{digit}{digit}{digit})
time                       ({digit}{digit}":"{digit}{digit}":"{digit}{digit})

versionNumber              ({digit}|".")*

comment                    {spaceNl}"(0"{space}
header                     {spaceNl}"(1"{space}
dimension                  {spaceNl}"(2"{space}
point                      {spaceNl}"(10"{space}
fluentFace                 {spaceNl}"(13"{space}
cell                       {spaceNl}"(12"{space}
zoneVariant1               {spaceNl}"(39"{space}
zoneVariant2               {spaceNl}"(45"{space}

unknownPeriodicFace        {spaceNl}"(17"{space}
periodicFace               {spaceNl}"(18"{space}
cellTree                   {spaceNl}"(58"{space}
faceTree                   {spaceNl}"(59"{space}
faceParents                {spaceNl}"(61"{space}

endOfSection               {space}")"{space}



 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%option stack

%x readComment
%x embeddedCommentState
%x readHeader
%x readDimension
%x readPoint
%x readPointHeader
%x readNumberOfPoints
%x readPointGroupData
%x readPointData
%x readPoints2D
%x readPoints3D
%x fluentFace
%x readFaceHeader
%x readNumberOfFaces
%x readFaceGroupData
%x readFaceData
%x readFacesMixed
%x readFacesUniform
%x cell
%x readCellHeader
%x readNumberOfCells
%x readCellGroupData
%x readCellData
%x readCellsMixed
%x readCellsUniform
%x zone
%x readZoneHeader
%x readZoneGroupData
%x readZoneData
%x readZoneBlock

%x periodicFace
%x cellTree
%x faceTree
%x faceParents

%x unknownBlock
%x embeddedUnknownBlock
%%

%{
    // Point data
    label pointGroupNumberOfComponents = 3;
    label pointI = 0; // index used for reading points

    // Face data
    label faceGroupElementType = -1;
    label faceI = 0;

    // Cell data
    label cellGroupElementType = -1;
    label cellI = 0;
%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

 /*                      ------ Reading control header ------                 */

{comment} {
        yy_push_state(readComment);
    }


<readComment>{quote}{text}{quote} {
    }


<readComment>{spaceNl}{endOfSection} {
        yy_pop_state();
    }

{header} {
        BEGIN(readHeader);
    }

<readHeader>{quote}{text}{quote} {
        Info << "Reading header: " << YYText() << endl;
    }


{dimension} {
        BEGIN(readDimension);
    }

<readDimension>{space}{label}{space} {
        IStringStream dimOfGridStream(YYText());

        dimensionOfGrid = readLabel(dimOfGridStream);

        Info << "Dimension of grid: " << dimensionOfGrid << endl;
    }


{point} {
        BEGIN(readPointHeader);
    }

<readPointHeader>{spaceNl}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfPoints);
    }

<readNumberOfPoints>{space}{hexLabel}{space}{labelList} {

        IStringStream numberOfPointsStream(YYText());

        nPoints = readHexLabel(numberOfPointsStream);

        Info << "Number of points: " << nPoints << endl;

        // set the size of the points list
        points.setSize(nPoints);

        // meaningless type skipped
        readLabel(numberOfPointsStream);

        // this dimension of grid may be checked against global dimension
        if (numberOfPointsStream)
        {
            // check dimension of grid
            readLabel(numberOfPointsStream);
        }
        else
        {
            Info << endl;
        }
    }

<readPointHeader>{spaceNl}{lbrac} {
        BEGIN(readPointGroupData);
    }

<readPointGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{labelList} {

        IStringStream pointGroupDataStream(YYText());

        // read point zone-ID, start and end-label
        // the indices will be used for checking later.
        pointGroupZoneID.append(readHexLabel(pointGroupDataStream));

        pointGroupStartIndex.append(readHexLabel(pointGroupDataStream));

        pointGroupEndIndex.append(readHexLabel(pointGroupDataStream));

        // point group type skipped
        readHexLabel(pointGroupDataStream);

        // In FOAM, indices start from zero - adjust
        pointI = pointGroupStartIndex.last() - 1;

        // reset number of components to default 
        pointGroupNumberOfComponents = 3;

        // read number of components in the vector
        if (pointGroupDataStream)
        {
            pointGroupNumberOfComponents = readLabel(pointGroupDataStream);
        }
    }

<readNumberOfPoints,readPointGroupData>{endOfSection} {
        BEGIN(readPointData);
    }

<readPointData>{spaceNl}{lbrac} {

        Info << "Reading points" << endl;

        if (pointGroupNumberOfComponents == 2)
        {
            yy_push_state(readPoints2D);
        }
        else
        {
            yy_push_state(readPoints3D);
        }
    }

<readPoints2D>{spaceNl}{scalarList} {

        IStringStream vertexXyzStream(YYText());

        // Note: coordinates must be read one at the time.  
        scalar x = readScalar(vertexXyzStream);
        scalar y = readScalar(vertexXyzStream);

        points[pointI] = point(x, y, 0);
        pointI++;
    }

<readPoints3D>{spaceNl}{scalarList} {

        IStringStream vertexXyzStream(YYText());

        // Note: coordinates must be read one at the time.  
        scalar x = readScalar(vertexXyzStream);
        scalar y = readScalar(vertexXyzStream);
        scalar z = readScalar(vertexXyzStream);

        points[pointI] = convertToMeters*point(x, y, z);
        pointI++;
    }

<readPoints2D,readPoints3D>{spaceNl}{endOfSection} {

        // check read of points
        if (pointI != pointGroupEndIndex.last())
        {
            Info << "problem with reading points: "
                << "start index: " << pointGroupStartIndex.last()
                << " end index: " << pointGroupEndIndex.last()
                << " last points read: " << pointI << endl;
        }

        yy_pop_state();
    }


{fluentFace} {
        BEGIN(readFaceHeader);
    }

<readFaceHeader>{spaceNl}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfFaces);
    }

<readNumberOfFaces>{space}{hexLabel}{space}{labelListElement}+ {

        IStringStream numberOfFacesStream(YYText());

        nFaces = readHexLabel(numberOfFacesStream);

        Info << "number of faces: " << nFaces << endl;

        faces.setSize(nFaces);
        owner.setSize(nFaces);
        neighbour.setSize(nFaces);

        // Meaningless type and element type not read
    }

<readFaceHeader>{spaceNl}{lbrac} {
        BEGIN(readFaceGroupData);
    }

<readFaceGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{hexLabelListElement}+ {

        IStringStream faceGroupDataStream(YYText());

        // read fluentFace zone-ID, start and end-label
        faceGroupZoneID.append(readHexLabel(faceGroupDataStream));

        // the indices will be used for checking later.
        faceGroupStartIndex.append(readHexLabel(faceGroupDataStream));

        faceGroupEndIndex.append(readHexLabel(faceGroupDataStream));

        // face group type
        readHexLabel(faceGroupDataStream);

        faceGroupElementType = readHexLabel(faceGroupDataStream);

        // In FOAM, indices start from zero - adjust
        faceI = faceGroupStartIndex.last() - 1;
    }

<readNumberOfFaces,readFaceGroupData>{spaceNl}{endOfSection} {
        BEGIN(readFaceData);
    }

<readFaceData>{spaceNl}{lbrac} {

        if (faceGroupElementType == 0)
        {
            Info << "Reading mixed faces" << endl;
            yy_push_state(readFacesMixed);
        }
        else
        {
            Info << "Reading uniform faces" << endl;
            yy_push_state(readFacesUniform);
        }
    }

<readFacesMixed>{spaceNl}{hexLabelList} {

        IStringStream mixedFaceStream(YYText());

        face& curFaceLabels = faces[faceI];

        // set size of label list
        curFaceLabels.setSize(readLabel(mixedFaceStream));

        forAll (curFaceLabels, i)
        {
            curFaceLabels[i] = readHexLabel(mixedFaceStream) - 1;
        }

        // read neighbour and owner. Neighbour comes first
        neighbour[faceI] = readHexLabel(mixedFaceStream) - 1;
        owner[faceI] = readHexLabel(mixedFaceStream) - 1;
        faceI++;
    }

<readFacesUniform>{spaceNl}{hexLabelList} {

        IStringStream mixedFaceStream(YYText());

        face& curFaceLabels = faces[faceI];

        // set size of label list. This is OK because in Fluent the type
        // for edge is 2, for triangle is 3 and for quad is 4
        curFaceLabels.setSize(faceGroupElementType);

        forAll (curFaceLabels, i)
        {
            curFaceLabels[i] = readHexLabel(mixedFaceStream) - 1;
        }

        // read neighbour and owner. Neighbour comes first
        neighbour[faceI] = readHexLabel(mixedFaceStream) - 1;
        owner[faceI] = readHexLabel(mixedFaceStream) - 1;
        faceI++;
    }

<readFacesMixed,readFacesUniform>{spaceNl}{endOfSection} {

        // check read of fluentFaces
        if (faceI != faceGroupEndIndex.last())
        {
            Info << "problem with reading fluentFaces: "
                << "start index: " << faceGroupStartIndex.last()
                << " end index: " << faceGroupEndIndex.last()
                << " last fluentFaces read: " << faceI << endl;
        }

        yy_pop_state();
    }


{cell} {
        BEGIN(readCellHeader);
    }

<readCellHeader>{spaceNl}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfCells);
    }

<readNumberOfCells>{space}{hexLabel}{space}{labelListElement}+ {

        IStringStream numberOfCellsStream(YYText());

        nCells = readHexLabel(numberOfCellsStream);

        Info << "Number of cells: " << nCells << endl;

        fluentCellModelID.setSize(nCells);

        // Meaningless type and element type not read
    }

<readCellHeader>{spaceNl}{lbrac} {
        BEGIN(readCellGroupData);
    }

<readCellGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel} {
        // Warning. This entry must be above the next one because of the lexing
        // rules. It is introduced to deal with the problem of reading
        // non-standard cell definition from Tgrid, which misses the type label.
        // 
        Info << "Tgrid syntax problem: " << YYText() << endl;
        IStringStream cellGroupDataStream(YYText());

        // read cell zone-ID, start and end-label
        cellGroupZoneID.append(readHexLabel(cellGroupDataStream));

        // the indices will be used for checking later.
        cellGroupStartIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupEndIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupType.append(readHexLabel(cellGroupDataStream));

        // Note. Potentially skip cell set if type is zero.
        // This entry does not exist in Tgrid files. 
        if (dimensionOfGrid == 2)
        {
            // Tgrid creating triangles
            cellGroupElementType = 1;
        }
        else
        {
            cellGroupElementType = 2;
        }

        // In FOAM, indices start from zero - adjust
        cellI = cellGroupStartIndex.last() - 1;

        if (cellGroupElementType != 0)
        {
            label lastIndex = cellGroupEndIndex.last();

            for (; cellI < lastIndex; cellI++)
            {
                fluentCellModelID[cellI] = cellGroupElementType;
            }
        }
    }

<readCellGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel} {
        // Warning. See above

        IStringStream cellGroupDataStream(YYText());

        // read cell zone-ID, start and end-label
        cellGroupZoneID.append(readHexLabel(cellGroupDataStream));

        // the indices will be used for checking later.
        cellGroupStartIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupEndIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupType.append(readHexLabel(cellGroupDataStream));

        // Note. Potentially skip cell set if type is zero.

        cellGroupElementType = readHexLabel(cellGroupDataStream);

        // In FOAM, indices start from zero - adjust
        cellI = cellGroupStartIndex.last() - 1;

        if (cellGroupElementType != 0)
        {
            Info << "Reading uniform cells" << endl;
            label lastIndex = cellGroupEndIndex.last();

            for (; cellI < lastIndex; cellI++)
            {
                fluentCellModelID[cellI] = cellGroupElementType;
            }
        }
    }

<readNumberOfCells,readCellGroupData>{endOfSection} {
        BEGIN(readCellData);
    }

<readCellData>{spaceNl}{lbrac} {
        Info << "Reading mixed cells" << endl;
        yy_push_state(readCellsMixed);
    }

<readCellsMixed>{spaceNl}{labelList} {

        IStringStream fluentCellModelIDStream(YYText());

        label cellID;
        while (fluentCellModelIDStream.read(cellID))
        {
            fluentCellModelID[cellI] = cellID;
            cellI++;
        }
    }

<readCellsMixed>{spaceNl}{endOfSection} {

        // check read of cells
        if (cellI != cellGroupEndIndex.last())
        {
            Info<< "problem with reading cells: "
                << "start index: " << cellGroupStartIndex.last()
                << " end index: " << cellGroupEndIndex.last()
                << " last cells read: " << cellI << endl;
        }

        yy_pop_state();
    }



{zoneVariant1} {
        BEGIN(readZoneHeader);
    }

{zoneVariant2} {
        BEGIN(readZoneHeader);
    }

<readZoneHeader>{spaceNl}{lbrac} {
        BEGIN(readZoneGroupData);
    }

<readZoneGroupData>{space}{label}{space}{word}{space}{word} {

        IStringStream zoneDataStream(YYText());

        // cell zone-ID not in hexadecimal!!! Inconsistency 
        label zoneID(readLabel(zoneDataStream));

        if (zoneID > maxZoneID - 1)
        {
            // resize the container
            maxZoneID = zoneID + zoneIDBuffer;

            patchTypeIDs.setSize(maxZoneID);
            patchNameIDs.setSize(maxZoneID);
        }

        zoneDataStream >> patchTypeIDs[zoneID];
        zoneDataStream >> patchNameIDs[zoneID];
    }

<readZoneGroupData>{space}{label}{space}{word}{space}{word}{space}{label} {
        // Fluent manual inconsistency, version 6.1.22
        IStringStream zoneDataStream(YYText());

        // cell zone-ID not in hexadecimal!!! Inconsistency 
        label zoneID(readLabel(zoneDataStream));

        if (zoneID > maxZoneID - 1)
        {
            // resize the container
            maxZoneID = zoneID + zoneIDBuffer;

            patchTypeIDs.setSize(maxZoneID);
            patchNameIDs.setSize(maxZoneID);
        }

        zoneDataStream >> patchTypeIDs[zoneID];
        zoneDataStream >> patchNameIDs[zoneID];
    }

<readZoneGroupData>{endOfSection} {
        BEGIN(readZoneData);
    }

<readZoneData>{spaceNl}{lbrac} {

        Info << "Reading zone data" << endl;
        yy_push_state(readZoneBlock);
    }

<readZoneBlock>{spaceNl}{schemeSymbolList} {
    }

<readZoneBlock>{lbrac} {
        Info << "Found unknown block in zone" << endl;
        yy_push_state(unknownBlock);
    }

<readZoneBlock>{endOfSection} {
        yy_pop_state();
    }



 /*             ------ Reading end of section and others ------               */

<readHeader,readDimension,readPointData,readFaceData,readCellData,readZoneData>{spaceNl}{endOfSection} {
        BEGIN(INITIAL);
    }

 /*    ------ Reading unknown type or non-standard comment ------             */


{lbrac}{label} {
        Info << "Found unknown block" << endl;
        yy_push_state(unknownBlock);
    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{schemeSymbol} {
    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{lbrac} {
        Info << "Embedded blocks in comment or unknown" << endl;
        yy_push_state(embeddedUnknownBlock);

    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{endOfSection} {
        Info << "Found end of section in unknown" << endl;
        yy_pop_state();
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{labelList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{hexLabelList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{scalarList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{schemeSymbolList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{text} {
    }


 /* ------ Ignore remaining space and \n s.  Any other characters are errors. */

<readPoints2D,readPoints3D>.|\n {
    Info << "ERROR! Do not understand characters: " << YYText() << endl;
    }
.|\n {}


 /*  ------ On EOF return to previous file, if none exists terminate. ------  */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include <fstream>

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("Fluent mesh file");
    argList::validOptions.insert("scale", "scale factor");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    scalar scaleFactor = 1.0;
    if (args.options().found("scale"))
    {
        scaleFactor = atof(args.options()["scale"].c_str());
    }

#   include "createTime.H"

    fileName fluentFile(args.args()[3]);
    std::ifstream fluentStream(fluentFile.c_str());

    if (!fluentStream)
    {
        FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
            << args.executable()
            << ": file " << fluentFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&fluentStream);
    while(lexer.yylex() != 0);

    Info << "\n\nFINISHED LEXING\n\n\n";
    // Lookup table giving number of vertices given a fluent cell type ID
    // Currently not used. 
    // label fluentModelNVertices[7] = {-1, 3, 4, 4, 8, 5, 6};

    // Lookup table giving number of vertices given a fluent cell type ID
    label fluentModelNFaces[7] = {-1, 3, 4, 4, 6, 5, 5};

    // Make a list of cell faces to be filled in for owner and neighbour

    labelListList cellFaces(nCells);

    labelList nFacesInCell(nCells, 0);

    forAll (cellFaces, cellI)
    {
        cellFaces[cellI].setSize(fluentModelNFaces[fluentCellModelID[cellI]]);
    }

    // fill in owner and neighbour

    forAll (owner, faceI)
    {
        if (owner[faceI] > -1)
        {
            label curCell = owner[faceI];
            cellFaces[curCell][nFacesInCell[curCell]] = faceI;

            nFacesInCell[curCell]++;
        }
    }

    forAll (neighbour, faceI)
    {
        if (neighbour[faceI] > -1)
        {
            label curCell = neighbour[faceI];
            cellFaces[curCell][nFacesInCell[curCell]] = faceI;

            nFacesInCell[curCell]++;
        }
    }

    // Construct shapes from face lists
    cellShapeList cellShapes(nCells);

    // Extrude mesh into 3-D

    Info << "dimension of grid: " << dimensionOfGrid << endl;
    faceList frontAndBackFaces;

    if (dimensionOfGrid == 2)
    {
        const scalar extrusionFactor = 0.01;

        boundBox box(max(points), min(points));

        const scalar zOffset = extrusionFactor*box.mag();

        // two-dimensional grid. Extrude in z-direction
        Info<< "Grid is 2-D. Extruding in z-direction by: "
            << 2*zOffset << endl;

        pointField oldPoints = points;

        const label pointOffset = oldPoints.size();

        points.setSize(2*pointOffset);

        label nNewPoints = 0;

        // Note: In order for the owner-neighbour rules to be right, the
        // points given by Fluent need to represent the FRONT plane of the
        // geometry. Therefore, the extrusion will be in -z direction
        // 
        forAll (oldPoints, pointI)
        {
            points[nNewPoints] = oldPoints[pointI];

            points[nNewPoints].z() = zOffset;

            nNewPoints++;
        }

        forAll (oldPoints, pointI)
        {
            points[nNewPoints] = oldPoints[pointI];

            points[nNewPoints].z() = -zOffset;

            nNewPoints++;
        }

        // 2-D shape recognition
        Info << "Creating shapes for 2-D cells"<< endl;

        // Set the number of empty faces
        frontAndBackFaces.setSize(2*nCells);

        forAll (fluentCellModelID, cellI)
        {
            switch (fluentCellModelID[cellI])
            {
                case 1:
                {
                    cellShapes[cellI] =
                        cellShape
                        (
                            extrudedTriangleCellShape
                            (
                                cellI,
                                cellFaces[cellI],
                                faces,
                                owner,
                                neighbour,
                                pointOffset,
                                frontAndBackFaces
                            )
                        );
                }
                break;

                case 3:
                {
                    cellShapes[cellI] =
                        cellShape
                        (
                            extrudedQuadCellShape
                            (
                                cellI,
                                cellFaces[cellI],
                                faces,
                                owner,
                                neighbour,
                                pointOffset,
                                frontAndBackFaces
                            )
                        );
                }
                break;

                default:
                {
                    FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
                        << "unrecognised  2-D cell shape: "
                        << fluentCellModelID[cellI]
                        << abort(FatalError);
                }
            }
        }

        // Create new faces
        forAll (faces, faceI)
        {

            if (faces[faceI].size() != 2)
            {
                FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
                    << "fluentMeshToFoam: a 2-D face defined with "
                    << faces[faceI].size() << " points." << endl;
            }

            labelList& newFace = faces[faceI];

            newFace.setSize(4);

            newFace[2] = newFace[1] + pointOffset;

            newFace[3] = newFace[0] + pointOffset;
        }

        // Create new cells from 2-D shapes
    }
    else
    {
        // 3-D shape recognition
        Info << "Creating shapes for 3-D cells"<< endl;
        forAll (fluentCellModelID, cellI)
        {
            if
            (
                fluentCellModelID[cellI] == 2    // tet
             || fluentCellModelID[cellI] == 4    // hex
             || fluentCellModelID[cellI] == 5    // pyramid
             || fluentCellModelID[cellI] == 6    // prism
            )
            {
                cellShapes[cellI] =
                    cellShape
                    (
                        create3DCellShape
                        (
                            cellI,
                            cellFaces[cellI],
                            faces,
                            owner,
                            neighbour,
                            fluentCellModelID[cellI]
                        )
                    );
            }
            else
            {
                FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
                    << "unrecognised 3-D cell shape: "
                    << fluentCellModelID[cellI]
                    << abort(FatalError);
            }
        }
    }

    // boundary faces are oriented such that the owner is zero and the face
    // area vector points into the domain. Turn them round before making patches

    forAll (faces, faceI)
    {
        if (owner[faceI] == -1)
        {
            // reverse face
            labelList oldFace = faces[faceI];

            forAllReverse(oldFace, i)
            {
                faces[faceI][oldFace.size() - i - 1] =
                    oldFace[i];
            }
        }
    }

    // Make boundary patches

    SLList<label>::iterator faceGroupZoneIDIter = faceGroupZoneID.begin();
    SLList<label>::iterator faceGroupStartIndexIter =
        faceGroupStartIndex.begin();
    SLList<label>::iterator faceGroupEndIndexIter = faceGroupEndIndex.begin();

    // Note. Not all groups of faces will be boundary patches.
    // Take care on construction
    faceListList boundary(faceGroupZoneID.size());
    wordList patchNames(faceGroupZoneID.size());
    wordList patchTypes(faceGroupZoneID.size());
    word defaultFacesType = emptyPolyPatch::typeName;
    wordList patchPhysicalTypes(faceGroupZoneID.size());

    label nPatches = 0;

    for
    (
        ;
        faceGroupZoneIDIter != faceGroupZoneID.end()
     && faceGroupStartIndexIter != faceGroupStartIndex.end()
     && faceGroupEndIndexIter != faceGroupEndIndex.end();
        ++faceGroupZoneIDIter,
        ++faceGroupStartIndexIter,
        ++faceGroupEndIndexIter
    )
    {
        // get face type and name
        const word& curPatchType = patchTypeIDs[faceGroupZoneIDIter()];

        const word& curPatchName = patchNameIDs[faceGroupZoneIDIter()];

        Info << "Creating patch for zone: " << faceGroupZoneIDIter()
            << " start: " << faceGroupStartIndexIter() 
            << " end: " << faceGroupEndIndexIter()
            << " type: " << curPatchType << " name: " << curPatchName << endl;

        // make patch labels
        label faceLabel = faceGroupStartIndexIter() - 1;

        faceList patchFaces(faceGroupEndIndexIter() - faceLabel);

        forAll (patchFaces, faceI)
        {
            if
            (
                faces[faceLabel].size() == 3
             || faces[faceLabel].size() == 4
            )
            {
                patchFaces[faceI] = face(faces[faceLabel]);
                faceLabel++;
            }
            else
            {
                FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
                    << "unrecognised face shape with "
                    << patchFaces[faceI].size() << " vertices"
                    << abort(FatalError);
            }
        }

        if 
        (
            curPatchType == "pressure"
         || curPatchType == "pressure-inlet"
         || curPatchType == "inlet-vent"
         || curPatchType == "intake-fan"
         || curPatchType == "pressure-outlet"
         || curPatchType == "exhaust-fan"
         || curPatchType == "outlet-vent"
         || curPatchType == "pressure-far-field"
         || curPatchType == "velocity-inlet"
         || curPatchType == "mass-flow-inlet"
         || curPatchType == "outflow"
        )
        {
            boundary[nPatches] = patchFaces;
            patchTypes[nPatches] = polyPatch::typeName;
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if (curPatchType == "wall")
        {
            boundary[nPatches] = patchFaces;
            patchTypes[nPatches] = wallPolyPatch::typeName;
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if
        (
            curPatchType == "symmetry"
         || curPatchType == "axis"
        )
        {
            boundary[nPatches] = patchFaces;
            patchTypes[nPatches] = symmetryPolyPatch::typeName;
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if
        (
            curPatchType == "interior"
         || curPatchType == "solid"
         || curPatchType == "fan"
         || curPatchType == "radiator"
         || curPatchType == "porous-jump"
        )
        {
            // do nothing - these are internal faces
            Info<< "Patch " << faceGroupZoneIDIter()
                << " contains solid or internal faces. Not added to boundary"
                << endl;
        }
        else if
        (
            curPatchType == ""
        )
        {
            Info<< "Patch " << faceGroupZoneIDIter()
                << ": Faces are defined but the "
                << "not created as a zone. Not added to boundary"
                << endl;
        }
        else
        {
            FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
                << "fluent patch type " << curPatchType << " not recognised."
                << abort(FatalError);
        }
            
    }

    if (dimensionOfGrid == 2)
    {
        Info << "Creating patch for front and back planes" << endl;

        boundary[nPatches] = frontAndBackFaces;
        patchTypes[nPatches] = emptyPolyPatch::typeName;
        patchNames[nPatches] = "frontAndBackPlanes";

        nPatches++;
    }

    // reset size of boundary
    boundary.setSize(nPatches);
    patchTypes.setSize(nPatches);
    patchNames.setSize(nPatches);

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::defaultRegion,
        patchNames,
        patchTypes,
        defaultFacesType,
        patchPhysicalTypes
    );

    // Scale the points
    points *= scaleFactor;

    // Mesh will auto-write on construction
    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        points,
        cellShapes,
        boundary,
        patchNames,
        patchTypes,
        defaultFacesType,
        patchPhysicalTypes
    );

    Info << "Checking mesh" << endl;
    pShapeMesh.checkMesh();

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing mesh" << endl;
    pShapeMesh.write();

    Info<< nl << "End" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of fluentMeshToFoam.L
 \* ------------------------------------------------------------------------- */
