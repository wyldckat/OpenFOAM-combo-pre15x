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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/
package FoamX.Options;

public class DocServerTab
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    private OptionsManager optionsMgr_;

    //--------------------------------------------------------------------------
    /** DocServerTab Constructor. */
    public DocServerTab(OptionsManager optionsMgr)
    {
        optionsMgr_ = optionsMgr;

        initComponents();

        if (optionsMgr_.containsKey("browser"))
        {
            browserField_.setText(optionsMgr_.getProperty("browser"));
        }
    }

    //--------------------------------------------------------------------------

    public void updateModel()
    {
        if (browserField_.getText().trim().length()> 0)
        {
            optionsMgr_.put("browser", browserField_.getText());
        }
        optionsMgr_ = null;
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        browserPanel_ = new javax.swing.JPanel();
        externalBrowser_ = new javax.swing.JCheckBox();
        browserLabel_ = new javax.swing.JLabel();
        browserField_ = new javax.swing.JTextField();
        serverPanel_ = new javax.swing.JPanel();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        browserPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        browserPanel_.setBorder(new javax.swing.border.TitledBorder(null, "Documentation Browser", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Dialog", 0, 10), java.awt.Color.black));
        externalBrowser_.setFont(new java.awt.Font("Dialog", 0, 10));
        externalBrowser_.setText("Use External Browser");
        externalBrowser_.setHorizontalTextPosition(javax.swing.SwingConstants.LEFT);
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(4, 4, 2, 4);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        browserPanel_.add(externalBrowser_, gridBagConstraints2);
        
        browserLabel_.setText("External Browser Executable");
        browserLabel_.setForeground(java.awt.Color.black);
        browserLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(2, 4, 4, 4);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        browserPanel_.add(browserLabel_, gridBagConstraints2);
        
        browserField_.setFont(new java.awt.Font("SansSerif", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(2, 4, 4, 4);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints2.weightx = 1.0;
        browserPanel_.add(browserField_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints1.weightx = 1.0;
        add(browserPanel_, gridBagConstraints1);
        
        serverPanel_.setBorder(new javax.swing.border.TitledBorder(null, "Documentation Server", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Dialog", 0, 10), java.awt.Color.black));
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(serverPanel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel browserPanel_;
    private javax.swing.JCheckBox externalBrowser_;
    private javax.swing.JLabel browserLabel_;
    private javax.swing.JTextField browserField_;
    private javax.swing.JPanel serverPanel_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
