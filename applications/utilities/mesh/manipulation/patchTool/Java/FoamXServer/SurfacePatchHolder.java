package FoamXServer;

/**
* FoamXServer/SurfacePatchHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/


// Surface patch description.
public final class SurfacePatchHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.SurfacePatch value = null;

  public SurfacePatchHolder ()
  {
  }

  public SurfacePatchHolder (FoamXServer.SurfacePatch initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.SurfacePatchHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.SurfacePatchHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.SurfacePatchHelper.type ();
  }

}
