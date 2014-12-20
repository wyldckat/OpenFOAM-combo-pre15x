package FoamXServer;


/**
* FoamXServer/FoamXAnyHelper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// -------------------------------------------------------------------------
abstract public class FoamXAnyHelper
{
  private static String  _id = "IDL:FoamXServer/FoamXAny:1.0";

  public static void insert (org.omg.CORBA.Any a, FoamXServer.FoamXAny that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static FoamXServer.FoamXAny extract (org.omg.CORBA.Any a)
  {
    return read (a.create_input_stream ());
  }

  private static org.omg.CORBA.TypeCode __typeCode = null;
  private static boolean __active = false;
  synchronized public static org.omg.CORBA.TypeCode type ()
  {
    if (__typeCode == null)
    {
      synchronized (org.omg.CORBA.TypeCode.class)
      {
        if (__typeCode == null)
        {
          if (__active)
          {
            return org.omg.CORBA.ORB.init().create_recursive_tc ( _id );
          }
          __active = true;
          org.omg.CORBA.StructMember[] _members0 = new org.omg.CORBA.StructMember [2];
          org.omg.CORBA.TypeCode _tcOf_members0 = null;
          _tcOf_members0 = FoamXServer.FoamXTypeHelper.type ();
          _members0[0] = new org.omg.CORBA.StructMember (
            "type",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_any);
          _members0[1] = new org.omg.CORBA.StructMember (
            "value",
            _tcOf_members0,
            null);
          __typeCode = org.omg.CORBA.ORB.init ().create_struct_tc (FoamXServer.FoamXAnyHelper.id (), "FoamXAny", _members0);
          __active = false;
        }
      }
    }
    return __typeCode;
  }

  public static String id ()
  {
    return _id;
  }

  public static FoamXServer.FoamXAny read (org.omg.CORBA.portable.InputStream istream)
  {
    FoamXServer.FoamXAny value = new FoamXServer.FoamXAny ();
    value.type = FoamXServer.FoamXTypeHelper.read (istream);
    value.value = istream.read_any ();
    return value;
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, FoamXServer.FoamXAny value)
  {
    FoamXServer.FoamXTypeHelper.write (ostream, value.type);
    ostream.write_any (value.value);
  }

}
