package FoamXServer;

/**
* FoamXServer/DimensionSetHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/

public final class DimensionSetHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.DimensionSet value = null;

  public DimensionSetHolder ()
  {
  }

  public DimensionSetHolder (FoamXServer.DimensionSet initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.DimensionSetHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.DimensionSetHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.DimensionSetHelper.type ();
  }

}
