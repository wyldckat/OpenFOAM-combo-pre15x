package FoamXServer.CaseServer;


/**
* FoamXServer/CaseServer/IGeometricFieldOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 02 April 2007 09:53:40 o'clock BST
*/


// ---------------------------------------------------------------------
public interface IGeometricFieldOperations 
{
  String name ();
  void getInternalFieldValue (FoamXServer.IDictionaryEntryHolder internalFieldValue) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
  void getPatchFieldParameters (String patchName, FoamXServer.IDictionaryEntryHolder patchFieldValue) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;

  // Returns true if this has been modified
  boolean modified ();
} // interface IGeometricFieldOperations