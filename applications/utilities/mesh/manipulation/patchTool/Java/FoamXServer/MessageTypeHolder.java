package FoamXServer;

/**
* FoamXServer/MessageTypeHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/

public final class MessageTypeHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.MessageType value = null;

  public MessageTypeHolder ()
  {
  }

  public MessageTypeHolder (FoamXServer.MessageType initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.MessageTypeHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.MessageTypeHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.MessageTypeHelper.type ();
  }

}
