package FoamXServer.CaseBrowser;

/**
* FoamXServer/CaseBrowser/ICaseBrowserHolder.java .
* Generated by the IDL-to-Java compiler (portable), version "3.1"
* from FoamX.idl
* 04 March 2005 01:10:48 o'clock GMT
*/


// ---------------------------------------------------------------------
public final class ICaseBrowserHolder implements org.omg.CORBA.portable.Streamable
{
  public FoamXServer.CaseBrowser.ICaseBrowser value = null;

  public ICaseBrowserHolder ()
  {
  }

  public ICaseBrowserHolder (FoamXServer.CaseBrowser.ICaseBrowser initialValue)
  {
    value = initialValue;
  }

  public void _read (org.omg.CORBA.portable.InputStream i)
  {
    value = FoamXServer.CaseBrowser.ICaseBrowserHelper.read (i);
  }

  public void _write (org.omg.CORBA.portable.OutputStream o)
  {
    FoamXServer.CaseBrowser.ICaseBrowserHelper.write (o, value);
  }

  public org.omg.CORBA.TypeCode _type ()
  {
    return FoamXServer.CaseBrowser.ICaseBrowserHelper.type ();
  }

}
