package FoamXServer;

/**
* FoamXServer/JobDescriptorHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// Job descriptor structure.
public final class JobDescriptorHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.JobDescriptor value = null;

  public JobDescriptorHolder ()
  {
  }

  public JobDescriptorHolder (FoamXServer.JobDescriptor initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.JobDescriptorHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.JobDescriptorHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.JobDescriptorHelper.type ();
  }

}
