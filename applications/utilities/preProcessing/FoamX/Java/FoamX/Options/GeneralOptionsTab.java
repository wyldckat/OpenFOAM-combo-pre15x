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
package FoamX.Options;

public class GeneralOptionsTab
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    private OptionsManager optionsMgr_;
    private OptionTableModel tableModel_;

    //--------------------------------------------------------------------------
    /** GeneralOptionsTab constructor. */
    public GeneralOptionsTab(OptionsManager optionsMgr)
    {
        optionsMgr_ = optionsMgr;
        tableModel_ = new OptionTableModel(optionsMgr);

        // Initialise the table model.
        //tableModel_.addParameter("Prompt for password at startup");

        initComponents();
    }

    //--------------------------------------------------------------------------

    public void updateModel()
    {
    }

    //--------------------------------------------------------------------------
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        panel_ = new javax.swing.JPanel();
        scrollPane_ = new javax.swing.JScrollPane();
        table_ = new javax.swing.JTable();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        panel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        panel_.setBorder(new javax.swing.border.EtchedBorder());
        table_.setModel(tableModel_);
        scrollPane_.setViewportView(table_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        panel_.add(scrollPane_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(10, 10, 10, 10);
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(panel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel panel_;
    private javax.swing.JScrollPane scrollPane_;
    private javax.swing.JTable table_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
