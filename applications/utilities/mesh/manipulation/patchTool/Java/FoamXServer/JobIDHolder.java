package FoamXServer;

/**
* FoamXServer/JobIDHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/


// Job ID.
public final class JobIDHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.JobID value = null;

  public JobIDHolder ()
  {
  }

  public JobIDHolder (FoamXServer.JobID initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.JobIDHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.JobIDHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.JobIDHelper.type ();
  }

}
