package FoamXServer;


/**
* FoamXServer/Point3Holder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/

public final class Point3Holder implements org.omg.CORBA.portable.Streamable
{
  public float value[] = null;

  public Point3Holder ()
  {
  }

  public Point3Holder (float[] initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.Point3Helper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.Point3Helper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.Point3Helper.type ();
  }

}
