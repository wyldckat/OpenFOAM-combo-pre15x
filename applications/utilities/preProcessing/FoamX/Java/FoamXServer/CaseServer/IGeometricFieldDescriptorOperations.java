package FoamXServer.CaseServer;


/**
* FoamXServer/CaseServer/IGeometricFieldDescriptorOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 02 April 2007 09:53:40 o'clock BST
*/


// ---------------------------------------------------------------------
public interface IGeometricFieldDescriptorOperations 
{
  FoamXServer.ITypeDescriptor typeDescriptor ();
  FoamXServer.ITypeDescriptor fieldTypeDescriptor ();
  FoamXServer.CaseServer.IGeometryDescriptor geometryDescriptor ();

  //- Field name (eg, "p").
  String name ();

  //- Field name (eg, "p").
  void name (String newName);

  //- Field description (eg, "Pressure").
  String description ();

  //- Field description (eg, "Pressure").
  void description (String newDescription);

  //- Field type name (eg, "scalar").
  String fieldTypeName ();

  //- Field type name (eg, "scalar").
  void fieldTypeName (String newFieldTypeName);

  //- Geometry type name (eg, "vol").
  String geometryTypeName ();

  //- Geometry type name (eg, "vol").
  void geometryTypeName (String newGeometryTypeName);

  //- Field dimensions.
  FoamXServer.DimensionSet dimensions ();

  //- Field dimensions.
  void dimensions (FoamXServer.DimensionSet newDimensions);
} // interface IGeometricFieldDescriptorOperations