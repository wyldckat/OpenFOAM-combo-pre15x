package FoamXServer.CaseServer;

/**
* FoamXServer/CaseServer/IApplicationClassHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// ---------------------------------------------------------------------
public final class IApplicationClassHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.CaseServer.IApplicationClass value = null;

  public IApplicationClassHolder ()
  {
  }

  public IApplicationClassHolder (FoamXServer.CaseServer.IApplicationClass initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.CaseServer.IApplicationClassHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.CaseServer.IApplicationClassHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.CaseServer.IApplicationClassHelper.type ();
  }

}
