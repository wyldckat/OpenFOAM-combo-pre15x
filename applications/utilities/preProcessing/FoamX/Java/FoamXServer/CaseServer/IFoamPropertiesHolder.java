package FoamXServer.CaseServer;

/**
* FoamXServer/CaseServer/IFoamPropertiesHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// ---------------------------------------------------------------------
public final class IFoamPropertiesHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.CaseServer.IFoamProperties value = null;

  public IFoamPropertiesHolder ()
  {
  }

  public IFoamPropertiesHolder (FoamXServer.CaseServer.IFoamProperties initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.CaseServer.IFoamPropertiesHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.CaseServer.IFoamPropertiesHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.CaseServer.IFoamPropertiesHelper.type ();
  }

}
