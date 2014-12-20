package FoamXServer.PatchToolServer;


/**
* FoamXServer/PatchToolServer/IPatchToolServerOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// ---------------------------------------------------------------------
public interface IPatchToolServerOperations 
{

  // Attributes.
  String caseRoot ();
  String caseName ();

  // Time-steps.
  String[] availableTimeSteps ();

  // Get current time
  String getTime () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Set time
  void setTime (String timeName, int timeIndex) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Read mesh and construct boundary mesh
  void read () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Read triangulated surface from file
  void readTriSurface (String fileName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Map onto mesh and write mesh
  void write () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Write triangulated surface file
  void writeTriSurface (String fileName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  //
  void getPoints (FoamXServer.FloatListHolder points) throws FoamXServer.FoamXError;
  int getNFaces () throws FoamXServer.FoamXError;
  int getNEdges () throws FoamXServer.FoamXError;
  void getEdges (FoamXServer.LongListHolder verts) throws FoamXServer.FoamXError;

  // Calculate feature edges
  void setFeatureEdges (double minCos) throws FoamXServer.FoamXError;

  // Compact points
  void getFeaturePoints (FoamXServer.FloatListHolder points) throws FoamXServer.FoamXError;

  // labels of featurePoints
  void getFeatureEdges (FoamXServer.LongListHolder edgeLabels) throws FoamXServer.FoamXError;

  // labels of edges that are features
  void getFeatureToEdge (FoamXServer.LongListHolder edgeLabels) throws FoamXServer.FoamXError;

  // From edge to featureEdge
  void getEdgeToFeature (FoamXServer.LongListHolder featLabels) throws FoamXServer.FoamXError;

  // Feature 'segments'. Arrays of indices into featureEdges
  void getFeatureSegments (FoamXServer.LongLongListHolder featureSegments) throws FoamXServer.FoamXError;

  // of an edge
  void setExtraEdges (int edgeI) throws FoamXServer.FoamXError;
  void getExtraEdges (FoamXServer.LongListHolder edgeLabels) throws FoamXServer.FoamXError;

  // Get patch of face
  int whichPatch (int faceI) throws FoamXServer.FoamXError;

  // Get index of patch by name
  int findPatchID (String patchName) throws FoamXServer.FoamXError;

  // Get all patches
  void getPatches (FoamXServer.SurfacePatchListHolder patches) throws FoamXServer.FoamXError;

  // Add patch
  void addPatch (String patchName) throws FoamXServer.FoamXError;

  // Delete patch
  void deletePatch (String patchName) throws FoamXServer.FoamXError;

  // Change patch type
  void changePatchType (String patchName, String patchType) throws FoamXServer.FoamXError;

  // Change patch faces
  void changeFaces (int[] patchIDs, FoamXServer.LongListHolder oldToNew) throws FoamXServer.FoamXError;

  // Return total number of triangles.
  int getNTris (int startFaceI, int nFaces, FoamXServer.LongListHolder nTris) throws FoamXServer.FoamXError;

  // Triangles returned as 3 consecutive vertices in tris.
  void triangulate (int startFaceI, int nFaces, int totalNTris, FoamXServer.LongListHolder triVerts) throws FoamXServer.FoamXError;

  // triangulation
  int getNPoints (int startFaceI, int nFaces) throws FoamXServer.FoamXError;

  // to global coords.
  void triangulateLocal (int startFaceI, int nFaces, int totalNTris, FoamXServer.LongListHolder triVerts, FoamXServer.LongListHolder localToGlobal) throws FoamXServer.FoamXError;

  // Flood filling
  void markFaces (int[] protectedEdges, int faceI, FoamXServer.BoolListHolder visited) throws FoamXServer.FoamXError;

  //- Lifetime management.
  void close ();
} // interface IPatchToolServerOperations