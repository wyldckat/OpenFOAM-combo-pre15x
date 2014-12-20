package FoamXServer;


/**
* FoamXServer/JobDescriptorHelper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/


// Job descriptor structure.
abstract public class JobDescriptorHelper
{
  private static String  _id = "IDL:FoamXServer/JobDescriptor:1.0";

  public static void insert (org.omg.CORBA.Any a, FoamXServer.JobDescriptor that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static FoamXServer.JobDescriptor extract (org.omg.CORBA.Any a)
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
          org.omg.CORBA.StructMember[] _members0 = new org.omg.CORBA.StructMember [19];
          org.omg.CORBA.TypeCode _tcOf_members0 = null;
          _tcOf_members0 = FoamXServer.JobIDHelper.type ();
          _members0[0] = new org.omg.CORBA.StructMember (
            "jobID",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_long);
          _members0[1] = new org.omg.CORBA.StructMember (
            "ppid",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_long);
          _members0[2] = new org.omg.CORBA.StructMember (
            "pgid",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[3] = new org.omg.CORBA.StructMember (
            "startDate",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[4] = new org.omg.CORBA.StructMember (
            "startTime",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[5] = new org.omg.CORBA.StructMember (
            "userName",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[6] = new org.omg.CORBA.StructMember (
            "foamVersion",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[7] = new org.omg.CORBA.StructMember (
            "code",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[8] = new org.omg.CORBA.StructMember (
            "argList",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[9] = new org.omg.CORBA.StructMember (
            "currentDir",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[10] = new org.omg.CORBA.StructMember (
            "rootDir",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[11] = new org.omg.CORBA.StructMember (
            "caseName",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_long);
          _members0[12] = new org.omg.CORBA.StructMember (
            "nProcs",
            _tcOf_members0,
            null);
          _tcOf_members0 = FoamXServer.JobIDHelper.type ();
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_sequence_tc (0, _tcOf_members0);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_alias_tc (FoamXServer.JobIDListHelper.id (), "JobIDList", _tcOf_members0);
          _members0[13] = new org.omg.CORBA.StructMember (
            "slaves",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_long);
          _members0[14] = new org.omg.CORBA.StructMember (
            "nCountedProcs",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().get_primitive_tc (org.omg.CORBA.TCKind.tk_double);
          _members0[15] = new org.omg.CORBA.StructMember (
            "cpuTime",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[16] = new org.omg.CORBA.StructMember (
            "endDate",
            _tcOf_members0,
            null);
          _tcOf_members0 = org.omg.CORBA.ORB.init ().create_string_tc (0);
          _members0[17] = new org.omg.CORBA.StructMember (
            "endTime",
            _tcOf_members0,
            null);
          _tcOf_members0 = FoamXServer.JobStatusHelper.type ();
          _members0[18] = new org.omg.CORBA.StructMember (
            "status",
            _tcOf_members0,
            null);
          __typeCode = org.omg.CORBA.ORB.init ().create_struct_tc (FoamXServer.JobDescriptorHelper.id (), "JobDescriptor", _members0);
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

  public static FoamXServer.JobDescriptor read (org.omg.CORBA.portable.InputStream istream)
  {
    FoamXServer.JobDescriptor value = new FoamXServer.JobDescriptor ();
    value.jobID = FoamXServer.JobIDHelper.read (istream);
    value.ppid = istream.read_long ();
    value.pgid = istream.read_long ();
    value.startDate = istream.read_string ();
    value.startTime = istream.read_string ();
    value.userName = istream.read_string ();
    value.foamVersion = istream.read_string ();
    value.code = istream.read_string ();
    value.argList = istream.read_string ();
    value.currentDir = istream.read_string ();
    value.rootDir = istream.read_string ();
    value.caseName = istream.read_string ();
    value.nProcs = istream.read_long ();
    value.slaves = FoamXServer.JobIDListHelper.read (istream);
    value.nCountedProcs = istream.read_long ();
    value.cpuTime = istream.read_double ();
    value.endDate = istream.read_string ();
    value.endTime = istream.read_string ();
    value.status = FoamXServer.JobStatusHelper.read (istream);
    return value;
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, FoamXServer.JobDescriptor value)
  {
    FoamXServer.JobIDHelper.write (ostream, value.jobID);
    ostream.write_long (value.ppid);
    ostream.write_long (value.pgid);
    ostream.write_string (value.startDate);
    ostream.write_string (value.startTime);
    ostream.write_string (value.userName);
    ostream.write_string (value.foamVersion);
    ostream.write_string (value.code);
    ostream.write_string (value.argList);
    ostream.write_string (value.currentDir);
    ostream.write_string (value.rootDir);
    ostream.write_string (value.caseName);
    ostream.write_long (value.nProcs);
    FoamXServer.JobIDListHelper.write (ostream, value.slaves);
    ostream.write_long (value.nCountedProcs);
    ostream.write_double (value.cpuTime);
    ostream.write_string (value.endDate);
    ostream.write_string (value.endTime);
    FoamXServer.JobStatusHelper.write (ostream, value.status);
  }

}
