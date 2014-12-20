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
import java.util.Hashtable;
import javax.swing.*;
import javax.swing.table.*;

import FoamX.App;
import FoamX.Editors.TypeEditor.TypeEditorPanel;
import FoamX.Exceptions.FoamXException;
import FoamX.Util.BusyCursor;
import FoamXServer.CaseServer.IFoamProperties;

public class ApplicationClassEditorDlg
    extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    private ApplicationClassModel appClassModel_;
    private Hashtable tabMap_;
    private ApplicationDetailsTab appDetailsTab_;
    private FieldsTab fieldsTab_;
    private TypeEditorPanel dictionariesTab_;
    private BoundaryTypesTab boundaryTypesTab_;
    private ModulesTab modulesTab_;
    private boolean cancelled_;
    private boolean readOnly_;

    //--------------------------------------------------------------------------
    /** Creates a new ApplicationClassEditor dialog based on the specified application class object. */
    public ApplicationClassEditorDlg(IFoamProperties foamSys, String appClassName, boolean readOnly) throws FoamXException
    {
        super(new javax.swing.JFrame(), true);    // Modal.

        try
        {
            // Set readonly flag.
            readOnly_ = readOnly;
            
            App.printMessage("Constructing Model...");
            // Create application class model.
            appClassModel_ = new ApplicationClassModel(foamSys);

            App.printMessage("Loading Model...");
            // Initialise the model object from the application class object.
            appClassModel_.loadApplicationClass(appClassName);

            // Create tabs.
            App.printMessage("Constructing App Details Tab...");
            appDetailsTab_    = new ApplicationDetailsTab(appClassModel_);
            App.printMessage("Constructing Field Tab...");
            fieldsTab_        = new FieldsTab(appClassModel_);
            App.printMessage("Constructing Boundary Types Tab...");
            boundaryTypesTab_ = new BoundaryTypesTab(appClassModel_);
            App.printMessage("Constructing Type Editor Tab...");
            dictionariesTab_  = new TypeEditorPanel(appClassModel_.getAppClassDescriptor());
            App.printMessage("Constructing Modules Tab...");
            modulesTab_       = new ModulesTab(appClassModel_);

            // Initialise components and add tabs.
            initComponents();
            tabbedPane.addTab("General", appDetailsTab_);
            tabbedPane.addTab("Fields", fieldsTab_);
            tabbedPane.addTab("Boundary Types", boundaryTypesTab_);
            tabbedPane.addTab("Dictionaries", dictionariesTab_);
            tabbedPane.addTab("Modules", modulesTab_);

            // Initialise tab map.
            tabMap_ = new Hashtable(5);
            tabMap_.put("General", appDetailsTab_);
            tabMap_.put("Fields", fieldsTab_);
            tabMap_.put("Boundary Types", boundaryTypesTab_);
            tabMap_.put("Dictionaries", dictionariesTab_);
            tabMap_.put("Modules", modulesTab_);

            cancelled_ = true;

            App.printMessage("Fin...");
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
            throw new FoamXException("Failed to construct ApplicationClassEditorDlg for class " + appClassName);
        }
    }

    //--------------------------------------------------------------------------

    public boolean wasCancelled()
    {
        return cancelled_;
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        tabbedPane = new javax.swing.JTabbedPane();
        buttonPanel = new javax.swing.JPanel();
        closeButton_ = new javax.swing.JButton();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setTitle("Application Class Editor");
        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setName("");
        setModal(true);
        setFont(new java.awt.Font("Dialog", 0, 10));
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                OnWindowClosing(evt);
            }
        });
        
        tabbedPane.setFont(new java.awt.Font("Dialog", 0, 10));
        tabbedPane.setAutoscrolls(true);
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(tabbedPane, gridBagConstraints1);
        
        buttonPanel.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        buttonPanel.setBorder(new javax.swing.border.EtchedBorder(javax.swing.border.EtchedBorder.RAISED));
        closeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        closeButton_.setLabel("Close");
        closeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnClose(evt);
            }
        });
        
        buttonPanel.add(closeButton_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.SOUTH;
        getContentPane().add(buttonPanel, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(650, 450));
        setLocation((screenSize.width-650)/2,(screenSize.height-450)/2);
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnClose(java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnClose

      closeWindow();

  }//GEN-LAST:event_OnClose

    //--------------------------------------------------------------------------

    /** Event handler for window closing event. */
  private void OnWindowClosing (java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_OnWindowClosing

      closeWindow();

  }//GEN-LAST:event_OnWindowClosing

    //--------------------------------------------------------------------------

    private void closeWindow()
    {
        try
        {
            if (!readOnly_)
            {
                int ret = JOptionPane.showConfirmDialog(this,
                                                         "Save changes to application class data?",
                                                         "FoamX...",
                                                         JOptionPane.YES_NO_CANCEL_OPTION,
                                                         JOptionPane.QUESTION_MESSAGE);
                if (ret == JOptionPane.YES_OPTION)
                {
                    if (!checkAndSaveData()) return;
                }
                else if (ret == JOptionPane.CANCEL_OPTION)
                {
                    return;
                }
            }

            setVisible(false);
            dispose();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private boolean checkAndSaveData()
    {
        // Show busy cursor.
        BusyCursor cursor = new BusyCursor(this);

        // Get each tab to update the model.
        appDetailsTab_.updateModel();
        fieldsTab_.updateModel();
        dictionariesTab_.updateModel();
        boundaryTypesTab_.updateModel();
        modulesTab_.updateModel();

        // Validate the model and see if it's happy.
        try
        {
            appClassModel_.validate();
        }
        catch (ApplicationClassExeption ex)
        {
            JOptionPane.showMessageDialog(this,
                                           "Error : " + ex.getMessage(),
                                           "FoamX...",
                                           JOptionPane.ERROR_MESSAGE);

            // Select offending tab.
            String tabName = ex.getTabName();
            if (tabMap_.containsKey(tabName))
            {
                Component tab = (Component)tabMap_.get(tabName);
                tabbedPane.setSelectedComponent(tab);
            }
            return false;
        }

        // Get model to serialise itself.
        appClassModel_.saveApplicationClass();

        // Indicate that we have committed the data.
        cancelled_ = false;

        return true;
    }

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTabbedPane tabbedPane;
    private javax.swing.JPanel buttonPanel;
    private javax.swing.JButton closeButton_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




