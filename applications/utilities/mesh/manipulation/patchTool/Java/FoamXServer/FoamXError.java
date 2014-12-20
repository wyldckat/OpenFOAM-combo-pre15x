package FoamXServer;


/**
* FoamXServer/FoamXError.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/

public final class FoamXError extends org.omg.CORBA.UserException
{
  public FoamXServer.ErrorCode errorCode = null;
  public String errorMessage = null;
  public String methodName = null;
  public String fileName = null;
  public int lineNo = (int)0;

  public FoamXError ()
  {
    super(FoamXErrorHelper.id());
  } // ctor

  public FoamXError (FoamXServer.ErrorCode _errorCode, String _errorMessage, String _methodName, String _fileName, int _lineNo)
  {
    super(FoamXErrorHelper.id());
    errorCode = _errorCode;
    errorMessage = _errorMessage;
    methodName = _methodName;
    fileName = _fileName;
    lineNo = _lineNo;
  } // ctor


  public FoamXError (String $reason, FoamXServer.ErrorCode _errorCode, String _errorMessage, String _methodName, String _fileName, int _lineNo)
  {
    super(FoamXErrorHelper.id() + "  " + $reason);
    errorCode = _errorCode;
    errorMessage = _errorMessage;
    methodName = _methodName;
    fileName = _fileName;
    lineNo = _lineNo;
  } // ctor

} // class FoamXError
