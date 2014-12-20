package FoamXServer;


/**
* FoamXServer/IDictionaryEntryOperations.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// -------------------------------------------------------------------------
public interface IDictionaryEntryOperations 
{

  // Reference to the type descriptor object for this entry.
  FoamXServer.ITypeDescriptor typeDescriptor ();

  // Current value for this (non compound) entry.
  FoamXServer.FoamXAny value ();

  // Current value for this (non compound) entry.
  void value (FoamXServer.FoamXAny newValue);
  void setValue (FoamXServer.FoamXAny value) throws FoamXServer.FoamXError;

  // Sub-elements for compound types.
  FoamXServer.IDictionaryEntry[] subElements ();
  int nSubElements () throws FoamXServer.FoamXError;
  boolean packedList () throws FoamXServer.FoamXError;

  // The current selection index
  int selection ();

  // The current selection index
  void selection (int newSelection);

  // entry object.
  void addElement (FoamXServer.IDictionaryEntryHolder subEntry) throws FoamXServer.FoamXError;

  // Remove element from list.
  void removeElement (FoamXServer.IDictionaryEntry subEntry) throws FoamXServer.FoamXError;

  // are found.
  void validate () throws FoamXServer.FoamXError, FoamXServer.ValidationError;

  // Returns true if this entry, or any sub-entries have been modified.
  boolean modified () throws FoamXServer.FoamXError;

  // Save method for root (dictionary) objects.
  void save () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError;
} // interface IDictionaryEntryOperations
