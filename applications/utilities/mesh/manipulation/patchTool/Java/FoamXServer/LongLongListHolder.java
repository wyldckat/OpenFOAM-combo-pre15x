package FoamXServer;


/**
* FoamXServer/LongLongListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 02 April 2007 10:12:21 o'clock BST
*/

public final class LongLongListHolder implements org.omg.CORBA.portable.Streamable
{
  public int value[][] = null;

  public LongLongListHolder ()
  {
  }

  public LongLongListHolder (int[][] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.LongLongListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.LongLongListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.LongLongListHelper.type ();
  }

}
