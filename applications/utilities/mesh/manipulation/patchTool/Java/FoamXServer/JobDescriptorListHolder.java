package FoamXServer;


/**
* FoamXServer/JobDescriptorListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/

public final class JobDescriptorListHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.JobDescriptor value[] = null;

  public JobDescriptorListHolder ()
  {
  }

  public JobDescriptorListHolder (FoamXServer.JobDescriptor[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.JobDescriptorListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.JobDescriptorListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.JobDescriptorListHelper.type ();
  }

}
