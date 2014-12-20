package FoamXServer;


/**
* FoamXServer/CaseDescriptorHelper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// Case descriptor structure.
abstract public class CaseDescriptorHelper
{
  private static String  _id = "IDL:FoamXServer/CaseDescriptor:1.0";

  public static void insert (org.omg.CORBA.Any a, FoamXServer.CaseDescriptor that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static FoamXServer.CaseDescriptor extract (org.omg.CORBA.Any a)
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
          org.omg.CORBA.StructMember[] _members0 = new org.omg.CORBA.StructMember [8];
          org.omg.CORBA.TypeCode _tcOf_members0 = null;
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[0] = new org.omg.CORBA.StructMember (
            "rootDir",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[1] = new org.omg.CORBA.StructMember (
            "rawRootDir",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[2] = new org.omg.CORBA.StructMember (
            "caseName",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[3] = new org.omg.CORBA.StructMember (
            "appClass",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_long);
          _members0[4] = new org.omg.CORBA.StructMember (
            "nProcs",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_boolean);
          _members0[5] = new org.omg.CORBA.StructMember (
            "managed",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_boolean);
          _members0[6] = new org.omg.CORBA.StructMember (
            "locked",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_boolean);
          _members0[7] = new org.omg.CORBA.StructMember (
            "error",
            _tcOf_members0,
            null);
          __typeCode = org.omg.CORBA.ORB.init ().create_struct_tc (FoamXServer.CaseDescriptorHelper.id (), "CaseDescriptor", _members0);
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

  public static FoamXServer.CaseDescriptor read (org.omg.CORBA.portable.InputStream istream)
  {
    FoamXServer.CaseDescriptor value = new FoamXServer.CaseDescriptor ();
    value.rootDir = istream.read_string ();
    value.rawRootDir = istream.read_string ();
    value.caseName = istream.read_string ();
    value.appClass = istream.read_string ();
    value.nProcs = istream.read_long ();
    value.managed = istream.read_boolean ();
    value.locked = istream.read_boolean ();
    value.error = istream.read_boolean ();
    return value;
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, FoamXServer.CaseDescriptor value)
  {
    ostream.write_string (value.rootDir);
    ostream.write_string (value.rawRootDir);
    ostream.write_string (value.caseName);
    ostream.write_string (value.appClass);
    ostream.write_long (value.nProcs);
    ostream.write_boolean (value.managed);
    ostream.write_boolean (value.locked);
    ostream.write_boolean (value.error);
  }

}
