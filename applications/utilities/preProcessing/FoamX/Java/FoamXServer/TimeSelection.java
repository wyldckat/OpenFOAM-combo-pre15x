package FoamXServer;


/**
* FoamXServer/TimeSelection.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* Thursday, March 13, 2003 4:23:04 PM GMT
*/

public class TimeSelection implements org.omg.CORBA.portable.IDLEntity
{
  private        int __value;
  private static int __size = 4;
  private static FoamXServer.TimeSelection[] __array = new FoamXServer.TimeSelection [__size];

  public static final int _TIME_LATESTTIME = 0;
  public static final FoamXServer.TimeSelection TIME_LATESTTIME = new FoamXServer.TimeSelection(_TIME_LATESTTIME);
  public static final int _TIME_FIRSTTIME = 1;
  public static final FoamXServer.TimeSelection TIME_FIRSTTIME = new FoamXServer.TimeSelection(_TIME_FIRSTTIME);
  public static final int _TIME_NOTIME = 2;
  public static final FoamXServer.TimeSelection TIME_NOTIME = new FoamXServer.TimeSelection(_TIME_NOTIME);
  public static final int _TIME_ALLTIMES = 3;
  public static final FoamXServer.TimeSelection TIME_ALLTIMES = new FoamXServer.TimeSelection(_TIME_ALLTIMES);

  public int value ()
  {
    return __value;
  }

  public static FoamXServer.TimeSelection from_int (int value)
  {
    if (value >= 0 && value < __size)
      return __array[value];
    else
      throw new org.omg.CORBA.BAD_PARAM ();
  }

  protected TimeSelection (int value)
  {
    __value = value;
    __array[__value] = this;
  }
} // class TimeSelection
