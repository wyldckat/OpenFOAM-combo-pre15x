package FoamXServer;


/**
* FoamXServer/JobID.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/


// Job ID.
public final class JobID implements org.omg.CORBA.portable.IDLEntity
{
  public String hostName = null;

  // Host name.
  public int processID = (int)0;

  public JobID ()
  {
  } // ctor

  public JobID (String _hostName, int _processID)
  {
    hostName = _hostName;
    processID = _processID;
  } // ctor

} // class JobID
