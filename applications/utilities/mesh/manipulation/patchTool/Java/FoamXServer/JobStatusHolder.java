package FoamXServer;

/**
* FoamXServer/JobStatusHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/

public final class JobStatusHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.JobStatus value = null;

  public JobStatusHolder ()
  {
  }

  public JobStatusHolder (FoamXServer.JobStatus initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.JobStatusHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.JobStatusHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.JobStatusHelper.type ();
  }

}