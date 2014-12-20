/*
** Support functions for writing binary FIELDVIEW unstructured data files.
*/

/* Include system stuff for I/O and string and exit functions. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Include the defines for the FV_* codes and wall_info values. */
#include "fv_reader_tags.h"


/* Don't change these - used by fv_encode_elem_header ! */
#define MAX_NUM_ELEM_FACES     6
#define BITS_PER_WALL  3
#define ELEM_TYPE_BIT_SHIFT    (MAX_NUM_ELEM_FACES*BITS_PER_WALL)


/*
** fv_encode_elem_header:  return an encoded binary element header
**
** Input:
**    elem_type:  integer element type as shown in fv_reader_tags.h
**    wall_info:  array of integer "wall" flags, one for each face of
**                the element.  The wall flags are used during streamline
**                calculation.  Currently, the only meaningful values are
**                A_WALL and NOT_A_WALL as shown in fv_reader_tags.h.
**                Streamlines are forced away from faces marked as
**                "A_WALL", by limiting velocity and position very near
**                the wall.
** Output:
**    Function return value is the encoded binary element header.
*/

#ifdef __STDC__
unsigned int fv_encode_elem_header (int elem_type, int wall_info[])
#else
unsigned int fv_encode_elem_header (elem_type, wall_info)
int elem_type;
int wall_info[];
#endif
{
    unsigned int header;
    int i, nfaces;

    switch (elem_type)
    {
        case FV_TET_ELEM_ID:
            header = (1 << ELEM_TYPE_BIT_SHIFT);
            nfaces = 4;
            break;
        case FV_HEX_ELEM_ID:
            header = (4 << ELEM_TYPE_BIT_SHIFT);
            nfaces = 6;
            break;
        case FV_PRISM_ELEM_ID:
            header = (3 << ELEM_TYPE_BIT_SHIFT);
            nfaces = 5;
            break;
        case FV_PYRA_ELEM_ID:
            header = (2 << ELEM_TYPE_BIT_SHIFT);
            nfaces = 5;
            break;
        default:
            fprintf(stderr, "ERROR:  Unknown element type\n");
            return 0;
    }

    for (i = 0; i < nfaces; i++)
    {
        unsigned int u = wall_info[i];
        if (u > A_WALL)
        {
            fprintf(stderr, "ERROR:  Bad wall value\n");
            return 0;
        }
        header |= (u << (i*BITS_PER_WALL));
    }
    return header;
}

/*
** fwrite_str80:  write out a string padded to 80 characters.
**
** Like fwrite, this returns the number of items written, which
** should be 80 if successful, and less than 80 if it failed.
*/
#ifdef __STDC__
size_t fwrite_str80 (char *str, FILE *fp)
#else
int fwrite_str80 (str, fp)
char *str;
FILE *fp;
#endif
{
    char cbuf[80];
    size_t len;
    int i;

    /* Most of this just to avoid garbage after the name. */
    len = strlen(str);
    strncpy(cbuf, str, len < 80 ? len : 80);

    for (i = len; i < 80; i++)
        cbuf[i] = '\0';  /* pad with zeros */

    return fwrite(cbuf, sizeof(char), 80, fp);
}


/*
** Sample program which writes out a single unstructured grid containing
** two hex elements, with pressure and velocity data at the nodes, and
** surface data for temperature and velocity on some of the boundaries.
** The FIELDVIEW unstructured file is written in binary format.
** For simplicity, no error checking is done on the calls to fwrite
** and fwrite_str80.
*/
#if 0	/***** CHANGE THIS TO "#if 1" TO RUN THE SAMPLE PROGRAM. *****/
int main()
{
    char *file_name = "dbl_hex.uns";
    FILE *outfp;
    int num_grids = 1;
    int num_face_types = 3;
    /*
    ** Note that one of the boundary type names is "wall."
    ** The boundary name "wall" has no special meaning in FIELDVIEW.
    ** Boundary types and element walls are independent pieces of
    ** information.  The only way to mark an element face as a wall
    ** (for streamline calculation) is with the wall_info array passed
    ** to fv_encode_elem_header.
    */
    static char *face_type_names[3] = { "bottom", "top", "wall" };
    /*
    ** Each boundary type is flagged with 1 or 0 depending on
    ** whether surface results are present or absent (see below).
    */
    static int results_flag[3]      =  {    1    ,   1  ,   0    };
    /*
    ** Each boundary type is flagged with 1 or 0 depending on
    ** whether surface normals can be calculated from a "right
    ** hand rule" (see below).
    */
    static int normals_flag[3]      =  {    1    ,   1  ,   0    };

    /*
    ** Note that vector variables are specified by a ';' and vector name
    ** following the first scalar name of 3 scalar components of the
    ** vector.  If writing 2-D results, the third component must still
    ** be provided here, and its values must be written in the variables
    ** section below (typically padded with zeros.)
    */
    int num_vars = 4;
    static char *var_names[4] = { "pressure", "uvel; velocity", "vvel", "wvel" };
    int num_bvars = 4;
    static char *bvar_names[4] = { "temperature", "uvel; velocity", "vvel", "wvel" };

    unsigned int elem_header;
    int grid, i;
    int ibuf[10];

    int nnodes = 12;	/* Number of nodes in the grid. */

    /* Constants. */
    static float consts[4] = { 1., 0., 0., 0. };

    /* XYZ coordinates of the nodes. */
    static float x[12] = { -1., -1., 1., 1., -1., -1., 1., 1., -1., -1., 1., 1. };
    static float y[12] = { -1., -1., -1., -1., 1., 1., 1., 1., 3., 3., 3., 3. };
    static float z[12] = { -1., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1. };

    /* hex1 and hex2 are hex elements, defined as an array of node numbers. */
    static int hex1[8] = {1,2,3,4,5,6,7,8};
    static int hex2[8] = {5,6,7,8,9,10,11,12};

    /*
    ** Face definitions for boundary faces.
    ** All faces have 4 vertices.  If the face is triangular,
    ** the last vertex should be zero.
    */
    static int bot_faces[4] = { 1,2,4,3 };
    static int top_faces[4] = { 9,10,12,11 };
    static int wall_faces[8][4] =
        { {1,2,6,5}, {5,6,10,9}, {3,4,8,7}, {7,8,12,11},
          {1,3,7,5}, {5,7,11,9}, {2,4,8,6}, {6,8,12,10} };

    /* Wall values for element faces. */
    static int hex1_walls[6] = { A_WALL, A_WALL, NOT_A_WALL,
	    NOT_A_WALL, A_WALL, A_WALL };
    static int hex2_walls[6] = { A_WALL, A_WALL, NOT_A_WALL,
	    NOT_A_WALL, A_WALL, A_WALL };

    /* 4 variables (pressure and velocity values) at the 12 grid nodes. */
    static float vars[4][12] =
      { {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11},
        {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2},
        {1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1},
        {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2} };

    /*
    ** 4 boundary variables (temperature and velocity values) defined on
    ** the single top face, and the the single bottom face.
    */
    static float top_bvars[4] = { 1.0, 2.0,4.0,2.5 };
    static float bot_bvars[4] = { 1.0, 4.5,3.0,3.0 };


    /* Open the file for binary write access. */
    if ((outfp = fopen(file_name, "wb")) == NULL)
    {
        perror ("Cannot open output file");
        exit(1);
    }

    /* Output the magic number. */
    ibuf[0] = FV_MAGIC;
    fwrite(ibuf, sizeof(int), 1, outfp);

    /* Output file header and version number. */
    fwrite_str80("FIELDVIEW", outfp);
    /*
    ** This version of the FIELDVIEW unstructured file is "2.5".
    ** This is written as two integers.
    */
    ibuf[0] = 2;
    ibuf[1] = 5;
    fwrite(ibuf, sizeof(int), 2, outfp);

    /* Output constants for time, fsmach, alpha and re. */
    fwrite(consts, sizeof(float), 4, outfp);

    /* Output the number of grids. */
    ibuf[0] = num_grids;
    fwrite(ibuf, sizeof(int), 1, outfp);

    /*
    ** Output the table of boundary types.
    ** Each boundary type is preceded by 2 integer flags.
    ** The first flag is an "surface results flag".
    ** A value of 1 means surface results will be present for this
    ** boundary type (if any boundary variables are specified in the
    ** boundary variable names table below).
    ** A value of 0 means no surface results will be present.
    ** The second flag indicates whether boundary faces of this type have
    ** consistent "clockness" for the purpose of calculating a surface
    ** normal.  A value of 1 means that all faces of this type are
    ** written following a "right hand rule" for clockness.  In other
    ** words, if the vertices are written on counter-clockwise:
    ** 4 --- 3
    ** |     |
    ** 1 --- 2
    ** then the normal to the face is pointing towards you (not away
    ** from you).  A value of 0 means that the faces do not have any
    ** consistent clockness.  The "clockness" of surface normals is
    ** only used for calculating certain special surface integrals
    ** that involve surface normals.  If the surface normals flag
    ** is 0, these special integrals will not be available.
    */
    ibuf[0] = num_face_types;
    fwrite(ibuf, sizeof(int), 1, outfp);
    for (i = 0; i < num_face_types; i++) {
	ibuf[0] = results_flag[i];
	ibuf[1] = normals_flag[i];
	fwrite(ibuf, sizeof(int), 2, outfp);
        fwrite_str80(face_type_names[i], outfp);
    }

    /* Output the table of variable names. */
    /* The number of variables can be zero. */
    ibuf[0] = num_vars;
    fwrite(ibuf, sizeof(int), 1, outfp);
    for (i = 0; i < num_vars; i++)
        fwrite_str80(var_names[i], outfp);

    /*
    ** Output the table of boundary variable names.
    ** Boundary variables are associated with boundary faces, rather than
    ** with grid nodes.
    ** FIELDVIEW will automatically append "[BNDRY]" to each name
    ** so boundary variables can be easily distinguished from ordinary
    ** (grid node) variables.
    ** The number of boundary variables can be different from the number
    ** of ordinary variables.  The number of boundary variables can be
    ** zero.
    */
    ibuf[0] = num_bvars;
    fwrite(ibuf, sizeof(int), 1, outfp);
    for (i = 0; i < num_vars; i++)
        fwrite_str80(bvar_names[i], outfp);

    /* Output grid data. */
    for (grid = 0; grid < num_grids; grid++)
    {
	/* Output the node definition section for this grid. */
        ibuf[0] = FV_NODES;
        ibuf[1] = nnodes;
        fwrite(ibuf, sizeof(int), 2, outfp);

	/*
	** Output the X, then Y, then Z node coordinates.
	** Note that all of the X coordinates are output before any of
	** the Y coordinates.
	*/
        fwrite(x, sizeof(float), nnodes, outfp);
        fwrite(y, sizeof(float), nnodes, outfp);
        fwrite(z, sizeof(float), nnodes, outfp);

	/*
        ** Output boundary faces of the 3 different types.
	** All faces have 4 vertices.  If the face is triangular,
	** the last vertex should be zero.
	** TIP: A single boundary type can be broken into several sections
	** if you prefer.  Also, boundary face sections do not have to
	** be in order.  You may have a section of 10 faces of type 3,
	** followed by a section of 20 faces of type 2, followed by a
	** section of 15 more faces of type 3.  Breaking a boundary
	** type into very many short sections is less efficient.  The
	** boundaries will require more memory and be somewhat
	** slower to calculate in FIELDVIEW.
	** 
	*/
        ibuf[0] = FV_FACES;
        ibuf[1] = 1;	/* first boundary type */
        ibuf[2] = 1;	/* number of faces of this type */
        fwrite(ibuf, sizeof(int), 3, outfp);
        fwrite(bot_faces, sizeof(int), 4, outfp);

        ibuf[0] = FV_FACES;
        ibuf[1] = 2;	/* second boundary type */
        ibuf[2] = 1;	/* number of faces of this type */
        fwrite(ibuf, sizeof(int), 3, outfp);
        fwrite(top_faces, sizeof(int), 4, outfp);

        ibuf[0] = FV_FACES;
        ibuf[1] = 3;	/* third boundary type */
        ibuf[2] = 8;	/* number of faces of this type */
        fwrite(ibuf, sizeof(int), 3, outfp);
        fwrite(wall_faces, sizeof(int), 8*4, outfp);
	/*
	** Start an elements section.
	** There may be as many elements sections as needed.
	** Each section may contain a single element type or a
	** mixture of element types.
	** For maximum efficiency, each section should contain
	** a significant percentage of the elements in the grid.
	** The most efficient case is a single section containing
	** all elements in the grid.
	*/
        ibuf[0] = FV_ELEMENTS;
        ibuf[1] = 0;  /* tet count */
        ibuf[2] = 2;  /* hex count */
        ibuf[3] = 0;  /* prism count */
        ibuf[4] = 0;  /* pyramid count */
        fwrite(ibuf, sizeof(int), 5, outfp);

        /* Write header for first element. */
        elem_header = fv_encode_elem_header(FV_HEX_ELEM_ID, hex1_walls);
	if (elem_header == 0)
	{
	    fprintf (stderr, "fv_encode_elem_header failed for first hex.\n");
	    exit(1);
	}
	fwrite (&elem_header, sizeof(elem_header), 1, outfp);

	/* Write node definition for first element. */
        fwrite(hex1, sizeof(int), 8, outfp);

        /* Write header for second element. */
        elem_header = fv_encode_elem_header(FV_HEX_ELEM_ID, hex2_walls);
	if (elem_header == 0)
	{
	    fprintf (stderr, "fv_encode_elem_header failed for second hex.\n");
	    exit(1);
	}
	fwrite (&elem_header, sizeof(elem_header), 1, outfp);

	/* Write node definition for second element. */
        fwrite(hex2, sizeof(int), 8, outfp);

        /*
	** Output the variables data.
	** You must write the section header even if the number
	** of variables is zero.
	*/
	ibuf[0] = FV_VARIABLES;
	fwrite(ibuf, sizeof(int), 1, outfp);
	/*
	** Note that all of the data for the first variable is output
	** before any of the data for the second variable.
	*/
	for (i = 0; i < num_vars; i++)
	    fwrite(vars[i], sizeof(float), nnodes, outfp);

        /*
	** Output the boundary variables data.
	** Remember that the Boundary Table above has a "surface results
	** flag" indicating which boundary types have surface results.
	** The data should be written in the same order as the faces in
	** the Boundary Faces section, skipping over faces whose boundary
	** type has a surface results flag of zero (false).
	** For each variable, you should write one number per boundary face.
	** You must write the section header even if the number of boundary
	** variables is zero.
	*/
	ibuf[0] = FV_BNDRY_VARS;
	fwrite(ibuf, sizeof(int), 1, outfp);
	/*
	** Note that all of the data for the first variable is output
	** before any of the data for the second variable.
	*/
	for (i = 0; i < num_bvars; i++) {
	    int num_faces;
	    /*
	    ** The data for the bottom face is written first for each
	    ** variable, because the bottom face was written first in the
	    ** "Boundary Faces" section.
	    ** The "wall" faces are skipped, because the surface results
	    ** flag for the wall boundary type was 0 (false) in the
	    ** Boundary Table section.
	    */
	    num_faces = 1;	/* number of bottom faces */
	    fwrite(&bot_bvars[i], sizeof(float), num_faces, outfp);
	    num_faces = 1;	/* number of top faces */
	    fwrite(&top_bvars[i], sizeof(float), num_faces, outfp);
	}
    }

    if (fclose(outfp) != 0)
    {
	perror ("Cannot close output file");
	exit(1);
    }

    return 0;
}
#endif	/* end commenting out the sample program */
