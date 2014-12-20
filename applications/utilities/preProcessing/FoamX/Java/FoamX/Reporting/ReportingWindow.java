/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/
package FoamX.Reporting;

import FoamX.App;

public class ReportingWindow
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** ReportingWindow constructor. */
    public ReportingWindow()
    {
        try
        {
            initComponents();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void printMessage(String message)
    {
        textArea_.append(message);
        textArea_.append("\n");

        // Force scrolling to end of text
        textArea_.setCaretPosition(textArea_.getText().length());
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        scrollPane_ = new javax.swing.JScrollPane();
        textArea_ = new javax.swing.JTextArea();
        
        setLayout(new java.awt.BorderLayout());
        
        textArea_.setEditable(false);
        textArea_.setFont(new java.awt.Font("Monospaced", 0, 10));
        scrollPane_.setViewportView(textArea_);
        
        add(scrollPane_, java.awt.BorderLayout.CENTER);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JScrollPane scrollPane_;
    private javax.swing.JTextArea textArea_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
