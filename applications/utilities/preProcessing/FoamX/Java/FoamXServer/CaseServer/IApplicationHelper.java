package FoamXServer.CaseServer;


/**
* FoamXServer/CaseServer/IApplicationHelper.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 02 April 2007 09:53:40 o'clock BST
*/


// ---------------------------------------------------------------------
abstract public class IApplicationHelper
{
  private static String  _id = "IDL:FoamXServer/CaseServer/IApplication:1.0";

  public static void insert (org.omg.CORBA.Any a, FoamXServer.CaseServer.IApplication that)
  {
    org.omg.CORBA.portable.OutputStream out = a.create_output_stream ();
    a.type (type ());
    write (out, that);
    a.read_value (out.create_input_stream (), type ());
  }

  public static FoamXServer.CaseServer.IApplication extract (org.omg.CORBA.Any a)
  {
    return read (a.create_input_stream ());
  }

  private static org.omg.CORBA.TypeCode __typeCode = null;
  synchronized public static org.omg.CORBA.TypeCode type ()
  {
    if (__typeCode == null)
    {
      __typeCode = org.omg.CORBA.ORB.init ().create_interface_tc (FoamXServer.CaseServer.IApplicationHelper.id (), "IApplication");
    }
    return __typeCode;
  }

  public static String id ()
  {
    return _id;
  }

  public static FoamXServer.CaseServer.IApplication read (org.omg.CORBA.portable.InputStream istream)
  {
    return narrow (istream.read_Object (_IApplicationStub.class));
  }

  public static void write (org.omg.CORBA.portable.OutputStream ostream, FoamXServer.CaseServer.IApplication value)
  {
    ostream.write_Object ((org.omg.CORBA.Object) value);
  }

  public static FoamXServer.CaseServer.IApplication narrow (org.omg.CORBA.Object obj)
  {
    if (obj == null)
      return null;
    else if (obj instanceof FoamXServer.CaseServer.IApplication)
      return (FoamXServer.CaseServer.IApplication)obj;
    else if (!obj._is_a (id ()))
      throw new org.omg.CORBA.BAD_PARAM ();
    else
    {
      org.omg.CORBA.portable.Delegate delegate = ((org.omg.CORBA.portable.ObjectImpl)obj)._get_delegate ();
      FoamXServer.CaseServer._IApplicationStub stub = new FoamXServer.CaseServer._IApplicationStub ();
      stub._set_delegate(delegate);
      return stub;
    }
  }

}