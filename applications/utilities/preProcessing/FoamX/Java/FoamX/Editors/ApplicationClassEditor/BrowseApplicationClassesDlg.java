/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/
package FoamX.Editors.ApplicationClassEditor;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import FoamX.App;
import FoamX.Util.BusyCursor;

import FoamXServer.ApplicationClassDescriptor;
import FoamXServer.CaseServer.IApplicationClassHolder;
import FoamXServer.CaseServer.IFoamProperties;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;

public class BrowseApplicationClassesDlg
    extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    private IFoamProperties foamSys_;
    private DefaultListModel listModel_;
    private boolean modified_;

    //--------------------------------------------------------------------------
    /** BrowseApplicationClassesDlg constructor. */
    public BrowseApplicationClassesDlg(IFoamProperties foamSys)
    {
        super(new javax.swing.JFrame(), true);

        try
        {
            modified_ = false;

            // Take a reference to the foam properties object.
            foamSys_ = foamSys;

            // Initialise the list of defined application classes.
            listModel_ = new DefaultListModel();
            initialiseAppClassList();

            // Initialise gui.
            initComponents();

            // Listen out for selection events.
            appClassList_.getSelectionModel().addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent evt)
                    {
                        OnSelectionChanged(evt);
                    }
                }
            );

            // Set the list renderer.
            appClassList_.setCellRenderer(new FoamX.Util.FoamXListRenderer());

            // Select first item.
            if (listModel_.getSize()> 0)
            {
                appClassList_.getSelectionModel().setLeadSelectionIndex(0);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public boolean getModified()
    {
        return modified_;
    }

    //--------------------------------------------------------------------------

    private void initialiseAppClassList()
    {
        try
        {
            listModel_.clear();
            ApplicationClassDescriptor[] appClasses = foamSys_.applicationClasses();

            // Add an entry for all defined application classes.
            for (int i = 0; i <appClasses.length; i++)
            {
                listModel_.addElement(new AppClassListItem(appClasses[i]));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void OnSelectionChanged(ListSelectionEvent evt)
    {
        try
        {
            // Get the node object for the current selection.
            if (!appClassList_.getSelectionModel().isSelectionEmpty())
            {
                // Get the selected application class.
                int              index = appClassList_.getSelectionModel().getLeadSelectionIndex();
                AppClassListItem item  = (AppClassListItem)listModel_.getElementAt(index);
                
                // Set the delete button status.
                delAppClassButton_.setEnabled(!item.isSystemClass());
        }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private String getNewAppClassName()
    {
        // Prompt user for a new Boundary Definition name.
        return (String)JOptionPane.showInputDialog(this,
                                                    "Enter New Application Class Name",
                                                    "FoamX...",
                                                    JOptionPane.PLAIN_MESSAGE, null, null, null);
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        appClassListPanel_ = new javax.swing.JPanel();
        scrollPane_ = new javax.swing.JScrollPane();
        appClassList_ = new javax.swing.JList();
        appClassButtonPanel_ = new javax.swing.JPanel();
        openAppClassButton_ = new javax.swing.JButton();
        newAppClassButton_ = new javax.swing.JButton();
        cloneAppClassButton_ = new javax.swing.JButton();
        delAppClassButton_ = new javax.swing.JButton();
        buttonPanel_ = new javax.swing.JPanel();
        buttonOK_ = new javax.swing.JButton();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setTitle("FoamX Application Classes");
        setName("appClassesDlg");
        setModal(true);
        setFont(new java.awt.Font("Dialog", 0, 10));
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                closeDialog(evt);
            }
        });
        
        appClassListPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        appClassListPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        appClassList_.setModel(listModel_);
        appClassList_.setFont(new java.awt.Font("Dialog", 0, 10));
        appClassList_.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                OnAppClassListMouseClicked(evt);
            }
        });
        
        scrollPane_.setViewportView(appClassList_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        appClassListPanel_.add(scrollPane_, gridBagConstraints2);
        
        appClassButtonPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints3;
        
        openAppClassButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        openAppClassButton_.setText("Open...");
        openAppClassButton_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnEditApplicationsClass(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 0;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints3.weightx = 1.0;
        appClassButtonPanel_.add(openAppClassButton_, gridBagConstraints3);
        
        newAppClassButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        newAppClassButton_.setText("New...");
        newAppClassButton_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnNewApplicationClass(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 1;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints3.weightx = 1.0;
        appClassButtonPanel_.add(newAppClassButton_, gridBagConstraints3);
        
        cloneAppClassButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        cloneAppClassButton_.setText("Clone...");
        cloneAppClassButton_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnCloneApplicationClass(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 2;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints3.weightx = 1.0;
        appClassButtonPanel_.add(cloneAppClassButton_, gridBagConstraints3);
        
        delAppClassButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        delAppClassButton_.setText("Delete...");
        delAppClassButton_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnDeleteApplicationClass(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 3;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints3.weightx = 1.0;
        appClassButtonPanel_.add(delAppClassButton_, gridBagConstraints3);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        appClassListPanel_.add(appClassButtonPanel_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(appClassListPanel_, gridBagConstraints1);
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        buttonOK_.setFont(new java.awt.Font("Dialog", 0, 10));
        buttonOK_.setText("Close");
        buttonOK_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnOK(evt);
            }
        });
        
        buttonPanel_.add(buttonOK_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.weightx = 1.0;
        getContentPane().add(buttonPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(380, 240));
        setLocation((screenSize.width-380)/2,(screenSize.height-240)/2);
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnCloneApplicationClass (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnCloneApplicationClass

        try
        {
            // Show busy cursor.
            BusyCursor cursor = new BusyCursor(this);

            // Clone the selected application class.
            if (!appClassList_.getSelectionModel().isSelectionEmpty())
            {
                // Get name of the selected application class.
                int              index = appClassList_.getSelectionModel().getLeadSelectionIndex();
                AppClassListItem item  = (AppClassListItem)listModel_.getElementAt(index);
                String           name  = item.getText();

                // Get new application class name.
                String appClassName = getNewAppClassName();
                if (appClassName != null && appClassName.trim().length()> 0)
                {
                    try
                    {
                        // Create a new application class object by cloning the selected one.
                        IApplicationClassHolder appClassHolder = new IApplicationClassHolder();
                        foamSys_.cloneApplicationClass(name, appClassName, appClassHolder);
                    }
                    catch (Exception ex)
                    {
                        Toolkit.getDefaultToolkit().beep();
                        JOptionPane.showMessageDialog(this,
                                          "Application class clone operation failed.",
                                          "FoamX...",
                                          JOptionPane.ERROR_MESSAGE);
                        return;
                    }

                    // Show the application class editor dialog for the new class.
                    ApplicationClassEditorDlg dlg = new ApplicationClassEditorDlg(foamSys_, appClassName, false);
                    dlg.show();

                    // If dialog was cancelled, tidy up.
                    if (dlg.wasCancelled())
                    {
                        // Delete application class.
                        foamSys_.deleteApplicationClass(appClassName);
                    }
                    else
                    {
                        // Add name to list.
                        AppClassListItem appClass = new AppClassListItem(appClassName);
                        listModel_.addElement(appClass);
                    }

                    // Set modified flag.
                    modified_ = true;

                    // Save user properties.
                    foamSys_.saveUserProperties();
                }
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnCloneApplicationClass

    //--------------------------------------------------------------------------

  private void OnAppClassListMouseClicked (java.awt.event.MouseEvent evt)
    {//GEN-FIRST:event_OnAppClassListMouseClicked

        try
        {
            // Show busy cursor.
            BusyCursor cursor = new BusyCursor(this);

            // Check for double click on a list item.
            if (evt.getClickCount() == 2)
            {
                if (!appClassList_.getSelectionModel().isSelectionEmpty())
                {
                    // Get name of application class to be edited.
                    int              index = appClassList_.getSelectionModel().getLeadSelectionIndex();
                    AppClassListItem item  = (AppClassListItem)listModel_.getElementAt(index);

                    // Show the application class editor dialog for this class.
                    ApplicationClassEditorDlg dlg = new ApplicationClassEditorDlg(foamSys_, item.getText(), item.isSystemClass());
                    dlg.show();

                    // Set modified flag.
                    modified_ = !dlg.wasCancelled();

                    // Save user properties.
                    foamSys_.saveUserProperties();
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnAppClassListMouseClicked

    //--------------------------------------------------------------------------

  private void OnDeleteApplicationClass (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnDeleteApplicationClass

        try
        {
            // Show busy cursor.
            BusyCursor cursor = new BusyCursor(this);

            // Remove the selected application classes.
            if (!appClassList_.getSelectionModel().isSelectionEmpty())
            {
                int minIndex = appClassList_.getSelectionModel().getMinSelectionIndex();
                int maxIndex = appClassList_.getSelectionModel().getMaxSelectionIndex();
                // Need to remove from the top so that it doesn't get its knickers in a twist.
                for (int i=maxIndex; i>= minIndex; i--)
                {
                    int              index = appClassList_.getSelectionModel().getLeadSelectionIndex();
                    AppClassListItem item  = (AppClassListItem)listModel_.getElementAt(index);

                    // Check whether this is a system class.
                    if (item.isSystemClass())
                    {
                        JOptionPane.showMessageDialog(this,
                                          "Can't delete system application class '" + item.getText() + ".",
                                          "FoamX...",
                                          JOptionPane.ERROR_MESSAGE);
                    }
                    else
                    {
                        if (JOptionPane.YES_OPTION == JOptionPane.showConfirmDialog(this,
                                                              "Are you sure you wish to delete user application class '" + item.getText() + "'?",
                                                              "FoamX...",
                                                              JOptionPane.YES_NO_OPTION,
                                                              JOptionPane.QUESTION_MESSAGE))
                        {
                            // Remove from kernel.
                            foamSys_.deleteApplicationClass(item.getText());

                            // Remove from list.
                            listModel_.removeElementAt(i);

                            // Set modified flag.
                            modified_ = true;
                        }
                    }
                }
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnDeleteApplicationClass

    //--------------------------------------------------------------------------

  private void OnNewApplicationClass (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnNewApplicationClass

        try
        {
            // Show busy cursor.
            BusyCursor cursor = new BusyCursor(this);

            String appClassName = getNewAppClassName();
            if (appClassName != null && appClassName.trim().length()> 0)
            {
                try
                {
                    // Create a new application class object.
                    IApplicationClassHolder appClassHolder = new IApplicationClassHolder();
                    foamSys_.addApplicationClass(appClassName, appClassHolder);
                }
                catch (Exception ex)
                {
                    Toolkit.getDefaultToolkit().beep();
                    JOptionPane.showMessageDialog(this,
                                      "Application class create operation failed.",
                                      "FoamX...",
                                      JOptionPane.ERROR_MESSAGE);
                    return;
                }

                // Show the application class editor dialog for this class.
                ApplicationClassEditorDlg dlg = new ApplicationClassEditorDlg(foamSys_, appClassName, false);
                dlg.show();

                // If dialog was cancelled, tidy up.
                if (dlg.wasCancelled())
                {
                    // Delete application class.
                    foamSys_.deleteApplicationClass(appClassName);
                }
                else
                {
                    // Add name to list.
                    AppClassListItem appClass = new AppClassListItem(appClassName);
                    listModel_.addElement(appClass);
                }

                // Set modified flag.
                modified_ = !dlg.wasCancelled();

                // Save user properties.
                foamSys_.saveUserProperties();
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnNewApplicationClass

    //--------------------------------------------------------------------------

  private void OnEditApplicationsClass (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnEditApplicationsClass

        try
        {
            // Show busy cursor.
            BusyCursor cursor = new BusyCursor(this);

            // Edit the selected application class.
            if (!appClassList_.getSelectionModel().isSelectionEmpty())
            {
                // Get name of application class to be edited.
                int              index = appClassList_.getSelectionModel().getLeadSelectionIndex();
                AppClassListItem item  = (AppClassListItem)listModel_.getElementAt(index);

                // Show the application class editor dialog for this class.
                ApplicationClassEditorDlg dlg = new ApplicationClassEditorDlg(foamSys_, item.getText(), item.isSystemClass());
                dlg.show();

                // Set modified flag.
                modified_ = !dlg.wasCancelled();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnEditApplicationsClass

    //--------------------------------------------------------------------------

  private void OnOK (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnOK

        try
        {
            if (modified_)
            {
                foamSys_.saveUserProperties();
            }

            setVisible(false);
            dispose();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

  }//GEN-LAST:event_OnOK

    //--------------------------------------------------------------------------
    /** Closes the dialog */
  private void closeDialog(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_closeDialog

        try
        {
            if (modified_)
            {
                foamSys_.saveUserProperties();
            }

            setVisible(false);
            dispose();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

  }//GEN-LAST:event_closeDialog

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel appClassListPanel_;
  private javax.swing.JScrollPane scrollPane_;
  private javax.swing.JList appClassList_;
  private javax.swing.JPanel appClassButtonPanel_;
  private javax.swing.JButton openAppClassButton_;
  private javax.swing.JButton newAppClassButton_;
  private javax.swing.JButton cloneAppClassButton_;
  private javax.swing.JButton delAppClassButton_;
  private javax.swing.JPanel buttonPanel_;
  private javax.swing.JButton buttonOK_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    static private class AppClassListItem
    implements FoamX.Util.FoamXListRenderer.FoamXListItem
    {
        //-----------------------------------------------------------------------------------------

        protected String appClassName_;
        protected boolean systemClass_;
        protected static Icon[] icons_;

        //--------------------------------------------------------------------------

        static
        {
            try
            {
                // Load icons.
                icons_ = new Icon[3];
                icons_[0] = App.getResources().getIcon("SystemAppClassImage");
                icons_[1] = App.getResources().getIcon("UserAppClassImage");
                icons_[2] = App.getResources().getIcon("OverridenAppClassImage");
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }

        //-----------------------------------------------------------------------------------------

        public AppClassListItem(ApplicationClassDescriptor desc)
        {
            appClassName_ = desc.name;
            systemClass_  = desc.systemClass;
        }

        //-----------------------------------------------------------------------------------------

        public AppClassListItem(String appClassName)
        {
            appClassName_ = appClassName;
            systemClass_  = false;
        }

        //-----------------------------------------------------------------------------------------

        public Icon getIcon()
        {
            return systemClass_ ? icons_[0] : icons_[1];
        }

        //-----------------------------------------------------------------------------------------

        public String getText()
        {
            return appClassName_;
        }

        //-----------------------------------------------------------------------------------------

        public boolean isSystemClass()
        {
            return systemClass_;
        }

        //-----------------------------------------------------------------------------------------
    }

    //--------------------------------------------------------------------------
}





