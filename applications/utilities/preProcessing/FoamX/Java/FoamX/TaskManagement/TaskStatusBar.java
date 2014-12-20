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
package FoamX.TaskManagement;

public class TaskStatusBar
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /** Creates new form TaskStatusBar */
    public TaskStatusBar()
    {
        initComponents();
    }

    //--------------------------------------------------------------------------

    String getStatusText()
    {
        return statusLabel_.getText();
    }
    void setStatusText(String text)
    {
        statusLabel_.setText(text);
    }

    int getProgress()
    {
        return progressBar_.getModel().getValue();
    }
    void setProgress(int value)
    {
        progressBar_.getModel().setValue(value);
    }

    int getProgressMin()
    {
        return progressBar_.getModel().getMinimum();
    }
    void   setProgressMin(int value)
    {
        progressBar_.getModel().setMinimum(value);
    }

    int getProgressMax()
    {
        return progressBar_.getModel().getMaximum();
    }
    void setProgressMax(int value)
    {
        progressBar_.getModel().setMaximum(value);
    }

    //--------------------------------------------------------------------------

    void showProgress(boolean visible)
    {
        progressBar_.setVisible(visible);
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        statusLabel_ = new javax.swing.JLabel();
        progressBar_ = new javax.swing.JProgressBar();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        statusLabel_.setText("Task");
        statusLabel_.setForeground(java.awt.Color.black);
        statusLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.insets = new java.awt.Insets(0, 10, 0, 5);
        add(statusLabel_, gridBagConstraints1);
        
        progressBar_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 1;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.insets = new java.awt.Insets(0, 5, 0, 10);
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints1.weightx = 1.0;
        add(progressBar_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel statusLabel_;
    private javax.swing.JProgressBar progressBar_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





