package FoamXServer;


/**
* FoamXServer/FoamXIOError.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 22 June 2005 14:28:29 o'clock BST
*/

public final class FoamXIOError extends org.omg.CORBA.UserException
{
  public String errorMessage = null;
  public String ioFileName = null;
  public int ioStartLineNo = (int)0;
  public int ioEndLineNo = (int)0;
  public String methodName = null;
  public String fileName = null;
  public int lineNo = (int)0;

  public FoamXIOError ()
  {
    super(FoamXIOErrorHelper.id());
  } // ctor

  public FoamXIOError (String _errorMessage, String _ioFileName, int _ioStartLineNo, int _ioEndLineNo, String _methodName, String _fileName, int _lineNo)
  {
    super(FoamXIOErrorHelper.id());
    errorMessage = _errorMessage;
    ioFileName = _ioFileName;
    ioStartLineNo = _ioStartLineNo;
    ioEndLineNo = _ioEndLineNo;
    methodName = _methodName;
    fileName = _fileName;
    lineNo = _lineNo;
  } // ctor


  public FoamXIOError (String $reason, String _errorMessage, String _ioFileName, int _ioStartLineNo, int _ioEndLineNo, String _methodName, String _fileName, int _lineNo)
  {
    super(FoamXIOErrorHelper.id() + "  " + $reason);
    errorMessage = _errorMessage;
    ioFileName = _ioFileName;
    ioStartLineNo = _ioStartLineNo;
    ioEndLineNo = _ioEndLineNo;
    methodName = _methodName;
    fileName = _fileName;
    lineNo = _lineNo;
  } // ctor

} // class FoamXIOError
