package FoamXServer.CasePostServer;

/**
* FoamXServer/CasePostServer/ICasePostServerHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// ---------------------------------------------------------------------
public final class ICasePostServerHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.CasePostServer.ICasePostServer value = null;

  public ICasePostServerHolder ()
  {
  }

  public ICasePostServerHolder (FoamXServer.CasePostServer.ICasePostServer initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.CasePostServer.ICasePostServerHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.CasePostServer.ICasePostServerHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.CasePostServer.ICasePostServerHelper.type ();
  }

}