package FoamXServer;


/**
* FoamXServer/ApplicationClassDescriptorListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/

public final class ApplicationClassDescriptorListHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.ApplicationClassDescriptor value[] = null;

  public ApplicationClassDescriptorListHolder ()
  {
  }

  public ApplicationClassDescriptorListHolder (FoamXServer.ApplicationClassDescriptor[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.ApplicationClassDescriptorListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.ApplicationClassDescriptorListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.ApplicationClassDescriptorListHelper.type ();
  }

}