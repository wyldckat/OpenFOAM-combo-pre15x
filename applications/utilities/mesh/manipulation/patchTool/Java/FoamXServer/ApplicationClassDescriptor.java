package FoamXServer;


/**
* FoamXServer/ApplicationClassDescriptor.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 02 April 2007 10:12:21 o'clock BST
*/


// Application class descriptor structure.
public final class ApplicationClassDescriptor implements org.omg.CORBA.portable.IDLEntity
{
  public String name = null;
  public String category = null;
  public boolean systemClass = false;

  public ApplicationClassDescriptor ()
  {
  } // ctor

  public ApplicationClassDescriptor (String _name, String _category, boolean _systemClass)
  {
    name = _name;
    category = _category;
    systemClass = _systemClass;
  } // ctor

} // class ApplicationClassDescriptor