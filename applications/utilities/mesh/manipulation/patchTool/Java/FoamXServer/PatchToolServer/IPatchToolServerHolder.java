package FoamXServer.PatchToolServer;

/**
* FoamXServer/PatchToolServer/IPatchToolServerHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// ---------------------------------------------------------------------
public final class IPatchToolServerHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.PatchToolServer.IPatchToolServer value = null;

  public IPatchToolServerHolder ()
  {
  }

  public IPatchToolServerHolder (FoamXServer.PatchToolServer.IPatchToolServer initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.PatchToolServer.IPatchToolServerHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.PatchToolServer.IPatchToolServerHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.PatchToolServer.IPatchToolServerHelper.type ();
  }

}
