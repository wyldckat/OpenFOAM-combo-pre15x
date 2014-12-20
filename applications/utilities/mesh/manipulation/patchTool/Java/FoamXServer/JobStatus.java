package FoamXServer;


/**
* FoamXServer/JobStatus.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from PatchTool.idl
* 03 March 2005 22:12:52 o'clock GMT
*/

public class JobStatus implements org.omg.CORBA.portable.IDLEntity
{
  private        int __value;
  private static int __size = 7;
  private static FoamXServer.JobStatus[] __array = new FoamXServer.JobStatus [__size];

  public static final int _JOB_UNDEFINED = 0;
  public static final FoamXServer.JobStatus JOB_UNDEFINED = new FoamXServer.JobStatus(_JOB_UNDEFINED);
  public static final int _JOB_LAUNCHING = 1;
  public static final FoamXServer.JobStatus JOB_LAUNCHING = new FoamXServer.JobStatus(_JOB_LAUNCHING);
  public static final int _JOB_RUNNING = 2;
  public static final FoamXServer.JobStatus JOB_RUNNING = new FoamXServer.JobStatus(_JOB_RUNNING);
  public static final int _JOB_STOPPING = 3;
  public static final FoamXServer.JobStatus JOB_STOPPING = new FoamXServer.JobStatus(_JOB_STOPPING);
  public static final int _JOB_SUSPENDED = 4;
  public static final FoamXServer.JobStatus JOB_SUSPENDED = new FoamXServer.JobStatus(_JOB_SUSPENDED);
  public static final int _JOB_FINISHED = 5;
  public static final FoamXServer.JobStatus JOB_FINISHED = new FoamXServer.JobStatus(_JOB_FINISHED);
  public static final int _JOB_ABORTED = 6;
  public static final FoamXServer.JobStatus JOB_ABORTED = new FoamXServer.JobStatus(_JOB_ABORTED);

  public int value ()
  {
    return __value;
  }

  public static FoamXServer.JobStatus from_int (int value)
  {
    if (value >= 0 && value < __size)
      return __array[value];
    else
      throw new org.omg.CORBA.BAD_PARAM ();
  }

  protected JobStatus (int value)
  {
    __value = value;
    __array[__value] = this;
  }
} // class JobStatus
