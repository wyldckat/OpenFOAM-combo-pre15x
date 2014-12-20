package FoamXServer;


/**
* FoamXServer/_IDictionaryEntryStub.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// -------------------------------------------------------------------------
public class _IDictionaryEntryStub extends org.omg.CORBA.portable.ObjectImpl implements FoamXServer.IDictionaryEntry
{


  // Reference to the type descriptor object for this entry.
  public FoamXServer.ITypeDescriptor typeDescriptor ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_typeDescriptor", true);
                $in = _invoke ($out);
                FoamXServer.ITypeDescriptor $result = FoamXServer.ITypeDescriptorHelper.read ($in);
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return typeDescriptor (        );
            } finally {
                _releaseReply ($in);
            }
  } // typeDescriptor


  // Current value for this (non compound) entry.
  public FoamXServer.FoamXAny value ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_value", true);
                $in = _invoke ($out);
                FoamXServer.FoamXAny $result = FoamXServer.FoamXAnyHelper.read ($in);
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return value (        );
            } finally {
                _releaseReply ($in);
            }
  } // value


  // Current value for this (non compound) entry.
  public void value (FoamXServer.FoamXAny newValue)
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_set_value", true);
                FoamXServer.FoamXAnyHelper.write ($out, newValue);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                value (newValue        );
            } finally {
                _releaseReply ($in);
            }
  } // value

  public void setValue (FoamXServer.FoamXAny value) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("setValue", true);
                FoamXServer.FoamXAnyHelper.write ($out, value);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                setValue (value        );
            } finally {
                _releaseReply ($in);
            }
  } // setValue


  // Sub-elements for compound types.
  public FoamXServer.IDictionaryEntry[] subElements ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_subElements", true);
                $in = _invoke ($out);
                FoamXServer.IDictionaryEntry $result[] = FoamXServer.DictionaryEntryListHelper.read ($in);
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return subElements (        );
            } finally {
                _releaseReply ($in);
            }
  } // subElements

  public int nSubElements () throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("nSubElements", true);
                $in = _invoke ($out);
                int $result = $in.read_long ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return nSubElements (        );
            } finally {
                _releaseReply ($in);
            }
  } // nSubElements

  public boolean packedList () throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("packedList", true);
                $in = _invoke ($out);
                boolean $result = $in.read_boolean ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return packedList (        );
            } finally {
                _releaseReply ($in);
            }
  } // packedList


  // The current selection index
  public int selection ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_selection", true);
                $in = _invoke ($out);
                int $result = $in.read_long ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return selection (        );
            } finally {
                _releaseReply ($in);
            }
  } // selection


  // The current selection index
  public void selection (int newSelection)
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_set_selection", true);
                $out.write_long (newSelection);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                selection (newSelection        );
            } finally {
                _releaseReply ($in);
            }
  } // selection


  // entry object.
  public void addElement (FoamXServer.IDictionaryEntryHolder subEntry) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("addElement", true);
                $in = _invoke ($out);
                subEntry.value = FoamXServer.IDictionaryEntryHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                addElement (subEntry        );
            } finally {
                _releaseReply ($in);
            }
  } // addElement


  // Remove element from list.
  public void removeElement (FoamXServer.IDictionaryEntry subEntry) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("removeElement", true);
                FoamXServer.IDictionaryEntryHelper.write ($out, subEntry);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                removeElement (subEntry        );
            } finally {
                _releaseReply ($in);
            }
  } // removeElement


  // are found.
  public void validate () throws FoamXServer.FoamXError, FoamXServer.ValidationError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("validate", true);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else if (_id.equals ("IDL:FoamXServer/ValidationError:1.0"))
                    throw FoamXServer.ValidationErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                validate (        );
            } finally {
                _releaseReply ($in);
            }
  } // validate


  // Returns true if this entry, or any sub-entries have been modified.
  public boolean modified () throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("modified", true);
                $in = _invoke ($out);
                boolean $result = $in.read_boolean ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return modified (        );
            } finally {
                _releaseReply ($in);
            }
  } // modified


  // Save method for root (dictionary) objects.
  public void save () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("save", true);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else if (_id.equals ("IDL:FoamXServer/FoamXIOError:1.0"))
                    throw FoamXServer.FoamXIOErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                save (        );
            } finally {
                _releaseReply ($in);
            }
  } // save

  // Type-specific CORBA::Object operations
  private static String[] __ids = {
    "IDL:FoamXServer/IDictionaryEntry:1.0"};

  public String[] _ids ()
  {
    return (String[])__ids.clone ();
  }

  private void readObject (java.io.ObjectInputStream s) throws java.io.IOException
  {
     String str = s.readUTF ();
     String[] args = null;
     java.util.Properties props = null;
     org.omg.CORBA.Object obj = org.omg.CORBA.ORB.init (args, props).string_to_object (str);
     org.omg.CORBA.portable.Delegate delegate = ((org.omg.CORBA.portable.ObjectImpl) obj)._get_delegate ();
     _set_delegate (delegate);
  }

  private void writeObject (java.io.ObjectOutputStream s) throws java.io.IOException
  {
     String[] args = null;
     java.util.Properties props = null;
     String str = org.omg.CORBA.ORB.init (args, props).object_to_string (this);
     s.writeUTF (str);
  }
} // class _IDictionaryEntryStub
