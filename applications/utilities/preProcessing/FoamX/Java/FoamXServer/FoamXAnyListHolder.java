package FoamXServer;


/**
* FoamXServer/FoamXAnyListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/

public final class FoamXAnyListHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.FoamXAny value[] = null;

  public FoamXAnyListHolder ()
  {
  }

  public FoamXAnyListHolder (FoamXServer.FoamXAny[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.FoamXAnyListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.FoamXAnyListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.FoamXAnyListHelper.type ();
  }

}
