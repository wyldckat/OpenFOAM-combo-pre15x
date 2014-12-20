package FoamXServer;


/**
* FoamXServer/FoamXType.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 22 June 2005 14:41:58 o'clock BST
*/


// -------------------------------------------------------------------------
public class FoamXType implements org.omg.CORBA.portable.IDLEntity
{
  private        int __value;
  private static int __size = 21;
  private static FoamXServer.FoamXType[] __array = new FoamXServer.FoamXType [__size];

  public static final int _Type_Undefined = 0;
  public static final FoamXServer.FoamXType Type_Undefined = new FoamXServer.FoamXType(_Type_Undefined);
  public static final int _Type_Boolean = 1;
  public static final FoamXServer.FoamXType Type_Boolean = new FoamXServer.FoamXType(_Type_Boolean);
  public static final int _Type_Label = 2;
  public static final FoamXServer.FoamXType Type_Label = new FoamXServer.FoamXType(_Type_Label);
  public static final int _Type_Scalar = 3;
  public static final FoamXServer.FoamXType Type_Scalar = new FoamXServer.FoamXType(_Type_Scalar);
  public static final int _Type_Char = 4;
  public static final FoamXServer.FoamXType Type_Char = new FoamXServer.FoamXType(_Type_Char);
  public static final int _Type_Word = 5;
  public static final FoamXServer.FoamXType Type_Word = new FoamXServer.FoamXType(_Type_Word);
  public static final int _Type_String = 6;
  public static final FoamXServer.FoamXType Type_String = new FoamXServer.FoamXType(_Type_String);
  public static final int _Type_RootDir = 7;
  public static final FoamXServer.FoamXType Type_RootDir = new FoamXServer.FoamXType(_Type_RootDir);
  public static final int _Type_RootAndCase = 8;
  public static final FoamXServer.FoamXType Type_RootAndCase = new FoamXServer.FoamXType(_Type_RootAndCase);
  public static final int _Type_CaseName = 9;
  public static final FoamXServer.FoamXType Type_CaseName = new FoamXServer.FoamXType(_Type_CaseName);
  public static final int _Type_HostName = 10;
  public static final FoamXServer.FoamXType Type_HostName = new FoamXServer.FoamXType(_Type_HostName);
  public static final int _Type_File = 11;
  public static final FoamXServer.FoamXType Type_File = new FoamXServer.FoamXType(_Type_File);
  public static final int _Type_Directory = 12;
  public static final FoamXServer.FoamXType Type_Directory = new FoamXServer.FoamXType(_Type_Directory);
  public static final int _Type_Time = 13;
  public static final FoamXServer.FoamXType Type_Time = new FoamXServer.FoamXType(_Type_Time);
  public static final int _Type_DimensionSet = 14;
  public static final FoamXServer.FoamXType Type_DimensionSet = new FoamXServer.FoamXType(_Type_DimensionSet);
  public static final int _Type_FixedList = 15;
  public static final FoamXServer.FoamXType Type_FixedList = new FoamXServer.FoamXType(_Type_FixedList);
  public static final int _Type_List = 16;
  public static final FoamXServer.FoamXType Type_List = new FoamXServer.FoamXType(_Type_List);
  public static final int _Type_Dictionary = 17;
  public static final FoamXServer.FoamXType Type_Dictionary = new FoamXServer.FoamXType(_Type_Dictionary);
  public static final int _Type_Selection = 18;
  public static final FoamXServer.FoamXType Type_Selection = new FoamXServer.FoamXType(_Type_Selection);
  public static final int _Type_Compound = 19;
  public static final FoamXServer.FoamXType Type_Compound = new FoamXServer.FoamXType(_Type_Compound);
  public static final int _Type_Field = 20;
  public static final FoamXServer.FoamXType Type_Field = new FoamXServer.FoamXType(_Type_Field);

  public int value ()
  {
    return __value;
  }

  public static FoamXServer.FoamXType from_int (int value)
  {
    if (value >= 0 && value < __size)
      return __array[value];
    else
      throw new org.omg.CORBA.BAD_PARAM ();
  }

  protected FoamXType (int value)
  {
    __value = value;
    __array[__value] = this;
  }
} // class FoamXType
