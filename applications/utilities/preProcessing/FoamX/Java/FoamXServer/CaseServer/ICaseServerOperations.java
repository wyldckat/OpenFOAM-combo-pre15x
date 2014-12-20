package FoamXServer.CaseServer;


/**
* FoamXServer/CaseServer/ICaseServerOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// ---------------------------------------------------------------------
public interface ICaseServerOperations 
{

  // Attributes.
  boolean managed ();

  // Attributes.
  void managed (boolean newManaged);
  String caseRoot ();
  String caseName ();

  // Application class for this case (read only).
  FoamXServer.CaseServer.IApplicationClass applicationClass ();

  // Foam system properties (read only).
  FoamXServer.CaseServer.IFoamProperties foamProperties ();

  // Time-steps.
  String[] availableTimeSteps ();

  // Get current time
  String getTime () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Set time
  void setTime (String timeName, int timeIndex) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Mesh.
  boolean meshDefined ();
  void readMesh () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void importMesh (String hostName, String rootDir, String caseName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Field values.
  void getFieldValues (String fieldName, FoamXServer.CaseServer.IGeometricFieldHolder fieldValues) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Patches.
  String[] patchNames ();
  void addPatch (String patchName, String boundaryType) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void deletePatch (String patchName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void deleteAllPatches () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Patch boundary conditions.
  void setPatchBoundaryType (String patchName, String boundaryType) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void getPatchBoundaryType (String patchName, org.omg.CORBA.StringHolder boundaryType) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Dictionaries.
  void getDictionary (String dictionaryName, boolean forceRead, FoamXServer.IDictionaryEntryHolder dictRoot) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Copy default dictionary into dictionaryName
  void copyDefaultDictionary (String dictionaryName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void readFile (String name, org.omg.CORBA.StringHolder contents) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void writeFile (String name, String contents) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Get modification date
  int fileModificationDate (String fileName) throws FoamXServer.FoamXIOError;

  //- Calculation control.
  int runCase (String arguments) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void killCase () throws FoamXServer.FoamXError;

  //- Validation.
  void validate () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError, FoamXServer.ValidationError;

  // Returns true if this has been modified
  boolean modified ();

  //- Persistence.
  void save () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError, FoamXServer.ValidationError;

  //- Lifetime management.
  void close ();
} // interface ICaseServerOperations