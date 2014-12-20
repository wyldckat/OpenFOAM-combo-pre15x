package FoamXServer.PatchToolServer;


/**
* FoamXServer/PatchToolServer/_IPatchToolServerStub.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 02 April 2007 10:12:21 o'clock BST
*/


// ---------------------------------------------------------------------
public class _IPatchToolServerStub extends org.omg.CORBA.portable.ObjectImpl implements FoamXServer.PatchToolServer.IPatchToolServer
{


  // Attributes.
  public String caseRoot ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_caseRoot", true);
                $in = _invoke ($out);
                String $result = $in.read_string ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return caseRoot (        );
            } finally {
                _releaseReply ($in);
            }
  } // caseRoot

  public String caseName ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_caseName", true);
                $in = _invoke ($out);
                String $result = $in.read_string ();
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return caseName (        );
            } finally {
                _releaseReply ($in);
            }
  } // caseName


  // Time-steps.
  public String[] availableTimeSteps ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("_get_availableTimeSteps", true);
                $in = _invoke ($out);
                String $result[] = FoamXServer.StringListHelper.read ($in);
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return availableTimeSteps (        );
            } finally {
                _releaseReply ($in);
            }
  } // availableTimeSteps


  // Get current time
  public String getTime () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getTime", true);
                $in = _invoke ($out);
                String $result = $in.read_string ();
                return $result;
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
                return getTime (        );
            } finally {
                _releaseReply ($in);
            }
  } // getTime


  // Set time
  public void setTime (String timeName, int timeIndex) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("setTime", true);
                $out.write_string (timeName);
                $out.write_long (timeIndex);
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
                setTime (timeName, timeIndex        );
            } finally {
                _releaseReply ($in);
            }
  } // setTime


  // Read mesh and construct boundary mesh
  public void read () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("read", true);
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
                read (        );
            } finally {
                _releaseReply ($in);
            }
  } // read


  // Read triangulated surface from file
  public void readTriSurface (String fileName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("readTriSurface", true);
                $out.write_string (fileName);
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
                readTriSurface (fileName        );
            } finally {
                _releaseReply ($in);
            }
  } // readTriSurface


  // Map onto mesh and write mesh
  public void write () throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("write", true);
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
                write (        );
            } finally {
                _releaseReply ($in);
            }
  } // write


  // Write triangulated surface file
  public void writeTriSurface (String fileName) throws FoamXServer.FoamXError, FoamXServer.FoamXIOError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("writeTriSurface", true);
                $out.write_string (fileName);
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
                writeTriSurface (fileName        );
            } finally {
                _releaseReply ($in);
            }
  } // writeTriSurface


  //
  public void getPoints (FoamXServer.FloatListHolder points) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getPoints", true);
                $in = _invoke ($out);
                points.value = FoamXServer.FloatListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getPoints (points        );
            } finally {
                _releaseReply ($in);
            }
  } // getPoints

  public int getNFaces () throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getNFaces", true);
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
                return getNFaces (        );
            } finally {
                _releaseReply ($in);
            }
  } // getNFaces

  public int getNEdges () throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getNEdges", true);
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
                return getNEdges (        );
            } finally {
                _releaseReply ($in);
            }
  } // getNEdges

  public void getEdges (FoamXServer.LongListHolder verts) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getEdges", true);
                $in = _invoke ($out);
                verts.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getEdges (verts        );
            } finally {
                _releaseReply ($in);
            }
  } // getEdges


  // Calculate feature edges
  public void setFeatureEdges (double minCos) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("setFeatureEdges", true);
                $out.write_double (minCos);
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
                setFeatureEdges (minCos        );
            } finally {
                _releaseReply ($in);
            }
  } // setFeatureEdges


  // Compact points
  public void getFeaturePoints (FoamXServer.FloatListHolder points) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getFeaturePoints", true);
                $in = _invoke ($out);
                points.value = FoamXServer.FloatListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getFeaturePoints (points        );
            } finally {
                _releaseReply ($in);
            }
  } // getFeaturePoints


  // labels of featurePoints
  public void getFeatureEdges (FoamXServer.LongListHolder edgeLabels) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getFeatureEdges", true);
                $in = _invoke ($out);
                edgeLabels.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getFeatureEdges (edgeLabels        );
            } finally {
                _releaseReply ($in);
            }
  } // getFeatureEdges


  // labels of edges that are features
  public void getFeatureToEdge (FoamXServer.LongListHolder edgeLabels) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getFeatureToEdge", true);
                $in = _invoke ($out);
                edgeLabels.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getFeatureToEdge (edgeLabels        );
            } finally {
                _releaseReply ($in);
            }
  } // getFeatureToEdge


  // From edge to featureEdge
  public void getEdgeToFeature (FoamXServer.LongListHolder featLabels) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getEdgeToFeature", true);
                $in = _invoke ($out);
                featLabels.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getEdgeToFeature (featLabels        );
            } finally {
                _releaseReply ($in);
            }
  } // getEdgeToFeature


  // Feature 'segments'. Arrays of indices into featureEdges
  public void getFeatureSegments (FoamXServer.LongLongListHolder featureSegments) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getFeatureSegments", true);
                $in = _invoke ($out);
                featureSegments.value = FoamXServer.LongLongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getFeatureSegments (featureSegments        );
            } finally {
                _releaseReply ($in);
            }
  } // getFeatureSegments


  // of an edge
  public void setExtraEdges (int edgeI) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("setExtraEdges", true);
                $out.write_long (edgeI);
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
                setExtraEdges (edgeI        );
            } finally {
                _releaseReply ($in);
            }
  } // setExtraEdges

  public void getExtraEdges (FoamXServer.LongListHolder edgeLabels) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getExtraEdges", true);
                $in = _invoke ($out);
                edgeLabels.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getExtraEdges (edgeLabels        );
            } finally {
                _releaseReply ($in);
            }
  } // getExtraEdges


  // Get patch of face
  public int whichPatch (int faceI) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("whichPatch", true);
                $out.write_long (faceI);
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
                return whichPatch (faceI        );
            } finally {
                _releaseReply ($in);
            }
  } // whichPatch


  // Get index of patch by name
  public int findPatchID (String patchName) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("findPatchID", true);
                $out.write_string (patchName);
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
                return findPatchID (patchName        );
            } finally {
                _releaseReply ($in);
            }
  } // findPatchID


  // Get all patches
  public void getPatches (FoamXServer.SurfacePatchListHolder patches) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getPatches", true);
                $in = _invoke ($out);
                patches.value = FoamXServer.SurfacePatchListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                getPatches (patches        );
            } finally {
                _releaseReply ($in);
            }
  } // getPatches


  // Add patch
  public void addPatch (String patchName) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("addPatch", true);
                $out.write_string (patchName);
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
                addPatch (patchName        );
            } finally {
                _releaseReply ($in);
            }
  } // addPatch


  // Delete patch
  public void deletePatch (String patchName) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("deletePatch", true);
                $out.write_string (patchName);
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
                deletePatch (patchName        );
            } finally {
                _releaseReply ($in);
            }
  } // deletePatch


  // Change patch type
  public void changePatchType (String patchName, String patchType) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("changePatchType", true);
                $out.write_string (patchName);
                $out.write_string (patchType);
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
                changePatchType (patchName, patchType        );
            } finally {
                _releaseReply ($in);
            }
  } // changePatchType


  // Change patch faces
  public void changeFaces (int[] patchIDs, FoamXServer.LongListHolder oldToNew) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("changeFaces", true);
                FoamXServer.LongListHelper.write ($out, patchIDs);
                $in = _invoke ($out);
                oldToNew.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                changeFaces (patchIDs, oldToNew        );
            } finally {
                _releaseReply ($in);
            }
  } // changeFaces


  // Return total number of triangles.
  public int getNTris (int startFaceI, int nFaces, FoamXServer.LongListHolder nTris) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getNTris", true);
                $out.write_long (startFaceI);
                $out.write_long (nFaces);
                $in = _invoke ($out);
                int $result = $in.read_long ();
                nTris.value = FoamXServer.LongListHelper.read ($in);
                return $result;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                return getNTris (startFaceI, nFaces, nTris        );
            } finally {
                _releaseReply ($in);
            }
  } // getNTris


  // Triangles returned as 3 consecutive vertices in tris.
  public void triangulate (int startFaceI, int nFaces, int totalNTris, FoamXServer.LongListHolder triVerts) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("triangulate", true);
                $out.write_long (startFaceI);
                $out.write_long (nFaces);
                $out.write_long (totalNTris);
                $in = _invoke ($out);
                triVerts.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                triangulate (startFaceI, nFaces, totalNTris, triVerts        );
            } finally {
                _releaseReply ($in);
            }
  } // triangulate


  // triangulation
  public int getNPoints (int startFaceI, int nFaces) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("getNPoints", true);
                $out.write_long (startFaceI);
                $out.write_long (nFaces);
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
                return getNPoints (startFaceI, nFaces        );
            } finally {
                _releaseReply ($in);
            }
  } // getNPoints


  // to global coords.
  public void triangulateLocal (int startFaceI, int nFaces, int totalNTris, FoamXServer.LongListHolder triVerts, FoamXServer.LongListHolder localToGlobal) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("triangulateLocal", true);
                $out.write_long (startFaceI);
                $out.write_long (nFaces);
                $out.write_long (totalNTris);
                $in = _invoke ($out);
                triVerts.value = FoamXServer.LongListHelper.read ($in);
                localToGlobal.value = FoamXServer.LongListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                triangulateLocal (startFaceI, nFaces, totalNTris, triVerts, localToGlobal        );
            } finally {
                _releaseReply ($in);
            }
  } // triangulateLocal


  // Flood filling
  public void markFaces (int[] protectedEdges, int faceI, FoamXServer.BoolListHolder visited) throws FoamXServer.FoamXError
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("markFaces", true);
                FoamXServer.LongListHelper.write ($out, protectedEdges);
                $out.write_long (faceI);
                $in = _invoke ($out);
                visited.value = FoamXServer.BoolListHelper.read ($in);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                if (_id.equals ("IDL:FoamXServer/FoamXError:1.0"))
                    throw FoamXServer.FoamXErrorHelper.read ($in);
                else
                    throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                markFaces (protectedEdges, faceI, visited        );
            } finally {
                _releaseReply ($in);
            }
  } // markFaces


  //- Lifetime management.
  public void close ()
  {
            org.omg.CORBA.portable.InputStream $in = null;
            try {
                org.omg.CORBA.portable.OutputStream $out = _request ("close", false);
                $in = _invoke ($out);
                return;
            } catch (org.omg.CORBA.portable.ApplicationException $ex) {
                $in = $ex.getInputStream ();
                String _id = $ex.getId ();
                throw new org.omg.CORBA.MARSHAL (_id);
            } catch (org.omg.CORBA.portable.RemarshalException $rm) {
                close (        );
            } finally {
                _releaseReply ($in);
            }
  } // close

  // Type-specific CORBA::Object operations
  private static String[] __ids = {
    "IDL:FoamXServer/PatchToolServer/IPatchToolServer:1.0"};

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
} // class _IPatchToolServerStub
