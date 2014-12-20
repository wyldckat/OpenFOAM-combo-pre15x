package FoamXServer.CaseBrowser;


/**
* FoamXServer/CaseBrowser/ICaseBrowserOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// ---------------------------------------------------------------------
public interface ICaseBrowserOperations 
{

  // Foam system properties (user data editable).
  FoamXServer.CaseServer.IFoamProperties foamProperties ();

  // Managed cases list.
  FoamXServer.CaseDescriptor[] cases ();

  // Running jobs list.
  FoamXServer.JobDescriptor[] runningJobs ();

  // Finished jobs list.
  FoamXServer.JobDescriptor[] finishedJobs ();

  // Get environment variable
  void getEnv (String envName, org.omg.CORBA.StringHolder hostName) throws FoamXServer.FoamXError;

  // Get the machines hostName
  void getHostName (org.omg.CORBA.StringHolder hostName) throws FoamXServer.FoamXError;

  // Get userName
  void getUserName (org.omg.CORBA.StringHolder userName) throws FoamXServer.FoamXError;

  // Get modification date
  int fileModificationDate (String fileName) throws FoamXServer.FoamXIOError;

  // Read file and store contents in string
  void readFile (String fileName, org.omg.CORBA.StringHolder contents) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Write contents of string as a file
  void writeFile (String fileName, String contents) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Utility execution.
  int invokeUtility (String hostName, String utilityName, String[] arguments, String logName, boolean backGround) throws FoamXServer.FoamXError;
  void refreshCaseList () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void addToCaseList (String rootDir) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Case manipulation.
  void openCase (FoamXServer.CaseDescriptor caseDesc) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void newCase (String rootDir, String caseName, String app) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void importCase (String rootDir, String caseName, String app) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void deleteCase (FoamXServer.CaseDescriptor caseDesc) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void cloneCase (FoamXServer.CaseDescriptor caseDesc, String newCaseRootDir, String newCaseName, String newAppClassName, String timeSel) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Try to resolve NameServer entry for caseServer
  boolean getCaseServerReference (String rootDir, String caseName, FoamXServer.CaseServer.ICaseServerHolder caseObj) throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;

  // start casePostServer
  void openCasePost (FoamXServer.CaseDescriptor caseDesc, int nProcs) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Try to resolve NameServer entry for casePostServer
  boolean getCasePostServerReference (String rootDir, String caseName, int nProcs, FoamXServer.CasePostServer.ICasePostServerHolder caseObj) throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;
  boolean caseLocked (FoamXServer.CaseDescriptor caseDesc) throws FoamXServer.FoamXError;
  void unlockCase (String rootDir, String caseName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void unlockCaseDescriptor (FoamXServer.CaseDescriptor caseDesc) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void addCase (String rootDir, String rawRootDir, String caseName, String app) throws FoamXServer.FoamXError;
  void caseOpen (String rootDir, String caseName) throws FoamXServer.FoamXError;
  boolean isCaseInError (FoamXServer.CaseDescriptor caseDesc) throws FoamXServer.FoamXError;
  void caseIsInError (FoamXServer.CaseDescriptor caseDesc) throws FoamXServer.FoamXError;

  // Process control
  void refreshJobsLists () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void checkRunningJobs () throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;
  void purgeRunningJobs () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void purgeFinishedJob (FoamXServer.JobID jobID) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void purgeFinishedJobs (int nDays) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Kill process.
  void kill (FoamXServer.JobID jobID) throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;

  // Suspend process
  void suspend (FoamXServer.JobID jobID) throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;

  // Continue process
  void cont (FoamXServer.JobID jobID) throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;

  //      now=false: next natural dump
  void end (FoamXServer.JobID jobID, String rootDir, String caseName, boolean now) throws FoamXServer.FoamXError, FoamXServer.FoamXSYSError;

  // Reset the job status
  void setStatus (FoamXServer.JobID jobID, FoamXServer.JobStatus jobStatus) throws FoamXServer.FoamXError;

  // Validation.
  void validate () throws FoamXServer.FoamXError, FoamXServer.ValidationError;

  // Persistence.
  void save () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError, FoamXServer.ValidationError;

  // Lifetime management.
  void close ();
} // interface ICaseBrowserOperations
