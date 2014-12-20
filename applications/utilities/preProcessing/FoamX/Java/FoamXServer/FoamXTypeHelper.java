package FoamXServer;


/**
* FoamXServer/FoamXTypeHelper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 02 April 2007 09:53:40 o'clock BST
*/


// -------------------------------------------------------------------------
abstract public class FoamXTypeHelper
{
  private static String  _id = "IDL:FoamXServer/FoamXType:1.0";

  public static void insert (org.omg.CORBA.Any a, FoamXServer.FoamXType that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static FoamXServer.FoamXType extract (org.omg.CORBA.Any a)
  {
    return read (a.create_input_stream ());
  }

  private static org.omg.CORBA.TypeCode __typeCode = null;
  synchronized public static org.omg.CORBA.TypeCode type ()
  {
    if (__typeCode == null)
    {
      __typeCode = org.omg.CORBA.ORB.init ().create_enum_tc (FoamXServer.FoamXTypeHelper.id (), "FoamXType", new String[] { "Type_Undefined", "Type_Boolean", "Type_Label", "Type_Scalar", "Type_Char", "Type_Word", "Type_String", "Type_RootDir", "Type_RootAndCase", "Type_CaseName", "Type_HostName", "Type_File", "Type_Directory", "Type_Time", "Type_DimensionSet", "Type_FixedList", "Type_List", "Type_Dictionary", "Type_Selection", "Type_Compound", "Type_Field"} );
    }
    return __typeCode;
  }

  public static String id ()
  {
    return _id;
  }

  public static FoamXServer.FoamXType read (org.omg.CORBA.portable.InputStream istream)
  {
    return FoamXServer.FoamXType.from_int (istream.read_long ());
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, FoamXServer.FoamXType value)
  {
    ostream.write_long (value.value ());
  }

}