package FoamXServer;

/**
* FoamXServer/CaseDescriptorHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// Case descriptor structure.
public final class CaseDescriptorHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.CaseDescriptor value = null;

  public CaseDescriptorHolder ()
  {
  }

  public CaseDescriptorHolder (FoamXServer.CaseDescriptor initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.CaseDescriptorHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.CaseDescriptorHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.CaseDescriptorHelper.type ();
  }

}