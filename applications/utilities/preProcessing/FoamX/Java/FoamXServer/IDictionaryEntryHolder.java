package FoamXServer;

/**
* FoamXServer/IDictionaryEntryHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 02 April 2007 09:53:40 o'clock BST
*/


// -------------------------------------------------------------------------
public final class IDictionaryEntryHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.IDictionaryEntry value = null;

  public IDictionaryEntryHolder ()
  {
  }

  public IDictionaryEntryHolder (FoamXServer.IDictionaryEntry initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.IDictionaryEntryHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.IDictionaryEntryHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.IDictionaryEntryHelper.type ();
  }

}