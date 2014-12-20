package FoamXServer;


/**
* FoamXServer/SurfacePatch.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// Surface patch description.
public final class SurfacePatch implements org.omg.CORBA.portable.IDLEntity
{
  public String name = null;

  // Name of patch
  public String geometricType = null;

  // Optional geometric type
  public int size = (int)0;

  // Size of patch
  public int start = (int)0;

  public SurfacePatch ()
  {
  } // ctor

  public SurfacePatch (String _name, String _geometricType, int _size, int _start)
  {
    name = _name;
    geometricType = _geometricType;
    size = _size;
    start = _start;
  } // ctor

} // class SurfacePatch