package FoamXServer;


/**
* FoamXServer/StringListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// -------------------------------------------------------------------------
public final class StringListHolder implements org.omg.CORBA.portable.Streamable
{
  public String value[] = null;

  public StringListHolder ()
  {
  }

  public StringListHolder (String[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.StringListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.StringListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.StringListHelper.type ();
  }

}