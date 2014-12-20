package FoamXServer;


/**
* FoamXServer/Point3Helper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/

abstract public class Point3Helper
{
  private static String  _id = "IDL:FoamXServer/Point3:1.0";

  public static void insert (org.omg.CORBA.Any a, float[] that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static float[] extract (org.omg.CORBA.Any a)
  {
    return read (a.create_input_stream ());
  }

  private static org.omg.CORBA.TypeCode __typeCode = null;
  synchronized public static org.omg.CORBA.TypeCode type ()
  {
    if (__typeCode == null)
    {
      __typeCode = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_float);
      __typeCode = org.omg.CORBA.ORB.init ().create_array_tc (3, __typeCode );
      __typeCode = org.omg.CORBA.ORB.init ().create_alias_tc (FoamXServer.Point3Helper.id (), "Point3", __typeCode);
    }
    return __typeCode;
  }

  public static String id ()
  {
    return _id;
  }

  public static float[] read (org.omg.CORBA.portable.InputStream istream)
  {
    float value[] = null;
    value = new float[3];
    for (int _o0 = 0;_o0 < (3); ++_o0)
    {
      value[_o0] = istream.read_float ();
    }
    return value;
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, float[] value)
  {
    if (value.length != (3))
      throw new org.omg.CORBA.MARSHAL (0, org.omg.CORBA.CompletionStatus.COMPLETED_MAYBE);
    for (int _i0 = 0;_i0 < (3); ++_i0)
    {
      ostream.write_float (value[_i0]);
    }
  }

}
