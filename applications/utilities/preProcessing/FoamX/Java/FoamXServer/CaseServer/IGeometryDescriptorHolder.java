package FoamXServer.CaseServer;

/**
* FoamXServer/CaseServer/IGeometryDescriptorHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// ---------------------------------------------------------------------
public final class IGeometryDescriptorHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.CaseServer.IGeometryDescriptor value = null;

  public IGeometryDescriptorHolder ()
  {
  }

  public IGeometryDescriptorHolder (FoamXServer.CaseServer.IGeometryDescriptor initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.CaseServer.IGeometryDescriptorHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.CaseServer.IGeometryDescriptorHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.CaseServer.IGeometryDescriptorHelper.type ();
  }

}
