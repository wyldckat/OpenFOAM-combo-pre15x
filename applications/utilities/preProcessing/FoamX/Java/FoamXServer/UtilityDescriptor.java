package FoamXServer;


/**
* FoamXServer/UtilityDescriptor.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// Utility descriptor structure.
public final class UtilityDescriptor implements org.omg.CORBA.portable.IDLEntity
{
  public String name = null;
  public String description = null;
  public String category = null;
  public boolean systemClass = false;

  // System defined application class.
  public boolean changesMesh = false;

  // Changes mesh
  public boolean changesFields = false;

  // uses; not any additional fields)
  public String clientBean = null;

  // Client side Java bean for this utility.
  public FoamXServer.ITypeDescriptor controlDict = null;

  public UtilityDescriptor ()
  {
  } // ctor

  public UtilityDescriptor (String _name, String _description, String _category, boolean _systemClass, boolean _changesMesh, boolean _changesFields, String _clientBean, FoamXServer.ITypeDescriptor _controlDict)
  {
    name = _name;
    description = _description;
    category = _category;
    systemClass = _systemClass;
    changesMesh = _changesMesh;
    changesFields = _changesFields;
    clientBean = _clientBean;
    controlDict = _controlDict;
  } // ctor

} // class UtilityDescriptor