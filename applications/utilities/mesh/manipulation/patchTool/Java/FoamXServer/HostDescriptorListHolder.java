package FoamXServer;


/**
* FoamXServer/HostDescriptorListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/

public final class HostDescriptorListHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.HostDescriptor value[] = null;

  public HostDescriptorListHolder ()
  {
  }

  public HostDescriptorListHolder (FoamXServer.HostDescriptor[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.HostDescriptorListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.HostDescriptorListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.HostDescriptorListHelper.type ();
  }

}
