package FoamXServer;


/**
* FoamXServer/DoubleListHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/

public final class DoubleListHolder implements org.omg.CORBA.portable.Streamable
{
  public double value[] = null;

  public DoubleListHolder ()
  {
  }

  public DoubleListHolder (double[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.DoubleListHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.DoubleListHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.DoubleListHelper.type ();
  }

}
