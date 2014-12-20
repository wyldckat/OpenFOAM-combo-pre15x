package FoamXServer.CaseServer;


/**
* FoamXServer/CaseServer/_IGeometricFieldStub.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 02 April 2007 09:53:40 o'clock BST
*/


// ---------------------------------------------------------------------
public class _IGeometricFieldStub extends org.omg.CORBA.portable.ObjectImpl implements FoamXServer.CaseServer.IGeometricField
{

  public String name ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_name", true);
                $in = _invoke ($out);
                String $result = $in.read_string ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return name (        );
            } finally {
                _releaseReply ($in);
            }
  } // name

  public void getInternalFieldValue (FoamXServer.IDictionaryEntryHolder internalFieldValue) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getInternalFieldValue", true);
                $in = _invoke ($out);
                internalFieldValue.value = FoamXServer.IDictionaryEntryHelper.read ($in);
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
                getInternalFieldValue (internalFieldValue        );
            } finally {
                _releaseReply ($in);
            }
  } // getInternalFieldValue

  public void getPatchFieldParameters (String patchName, FoamXServer.IDictionaryEntryHolder patchFieldValue) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getPatchFieldParameters", true);
                $out.write_string (patchName);
                $in = _invoke ($out);
                patchFieldValue.value = FoamXServer.IDictionaryEntryHelper.read ($in);
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
                getPatchFieldParameters (patchName, patchFieldValue        );
            } finally {
                _releaseReply ($in);
            }
  } // getPatchFieldParameters


  // Returns true if this has been modified
  public boolean modified ()
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
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return modified (        );
            } finally {
                _releaseReply ($in);
            }
  } // modified

  // Type-specific CORBA::Object operations
  private static String[] __ids = {
    "IDL:FoamXServer/CaseServer/IGeometricField:1.0"};

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
} // class _IGeometricFieldStub