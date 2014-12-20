package FoamXServer;


/**
* FoamXServer/JobIDListHelper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 02 April 2007 10:12:21 o'clock BST
*/

abstract public class JobIDListHelper
{
  private static String  _id = "IDL:FoamXServer/JobIDList:1.0";

  public static void insert (org.omg.CORBA.Any a, FoamXServer.JobID[] that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static FoamXServer.JobID[] extract (org.omg.CORBA.Any a)
  {
    return read (a.create_input_stream ());
  }

  private static org.omg.CORBA.TypeCode __typeCode = null;
  synchronized public static org.omg.CORBA.TypeCode type ()
  {
    if (__typeCode == null)
    {
      __typeCode = FoamXServer.JobIDHelper.type ();
      __typeCode = org.omg.CORBA.ORB.init ().create_sequence_tc (0, __typeCode);
      __typeCode = org.omg.CORBA.ORB.init ().create_alias_tc (FoamXServer.JobIDListHelper.id (), "JobIDList", __typeCode);
    }
    return __typeCode;
  }

  public static String id ()
  {
    return _id;
  }

  public static FoamXServer.JobID[] read (org.omg.CORBA.portable.InputStream istream)
  {
    FoamXServer.JobID value[] = null;
    int _len0 = istream.read_long ();
    value = new FoamXServer.JobID[_len0];
    for (int _o1 = 0;_o1 < value.length; ++_o1)
      value[_o1] = FoamXServer.JobIDHelper.read (istream);
    return value;
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, FoamXServer.JobID[] value)
  {
    ostream.write_long (value.length);
    for (int _i0 = 0;_i0 < value.length; ++_i0)
      FoamXServer.JobIDHelper.write (ostream, value[_i0]);
  }

}
