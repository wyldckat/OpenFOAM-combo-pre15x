package FoamXServer;

/**
* FoamXServer/StringPairHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/

public final class StringPairHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.StringPair value = null;

  public StringPairHolder ()
  {
  }

  public StringPairHolder (FoamXServer.StringPair initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.StringPairHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.StringPairHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.StringPairHelper.type ();
  }

}
