package FoamXServer.CaseServer;


/**
* FoamXServer/CaseServer/IApplicationClassOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// ---------------------------------------------------------------------
public interface IApplicationClassOperations 
{

  //- Short name for this application class.
  String name ();

  //- Short name for this application class.
  void name (String newName);

  //- Application class description.
  String description ();

  //- Application class description.
  void description (String newDescription);

  //- Application class category.
  String category ();

  //- Application class category.
  void category (String newCategory);

  //- List of modules required to pre-process this application class.
  String[] modules ();

  //- List of modules required to pre-process this application class.
  void modules (String[] newModules);

  //- System or user defined application class.
  boolean systemClass ();

  //- List of defined fields.
  String[] fields ();
  void getField (String fieldName, FoamXServer.CaseServer.IGeometricFieldDescriptorHolder fieldDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void findField (String fieldName, FoamXServer.CaseServer.IGeometricFieldDescriptorHolder fieldDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void addField (String fieldName, FoamXServer.CaseServer.IGeometricFieldDescriptorHolder fieldDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void deleteField (String fieldName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Boundary types.
  String[] boundaryTypes ();

  // List of defined boundary types.
  void getBoundaryType (String boundaryTypeName, FoamXServer.CaseServer.IBoundaryTypeDescriptorHolder boundaryDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void findBoundaryType (String boundaryTypeName, FoamXServer.CaseServer.IBoundaryTypeDescriptorHolder boundaryDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void addBoundaryType (String boundaryTypeName, FoamXServer.CaseServer.IBoundaryTypeDescriptorHolder boundaryDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void deleteBoundaryType (String boundaryTypeName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Dictionaries.
  String[] dictionaries ();

  // List of defined dictionaries.
  void getDictionary (String dictName, FoamXServer.ITypeDescriptorHolder dictTypeDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void addDictionary (String dictName, FoamXServer.ITypeDescriptorHolder dictTypeDescriptor) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void deleteDictionary (String dictName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Validation.
  void validate () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError, FoamXServer.ValidationError;

  // Persistence.
  void save () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError, FoamXServer.ValidationError;
} // interface IApplicationClassOperations