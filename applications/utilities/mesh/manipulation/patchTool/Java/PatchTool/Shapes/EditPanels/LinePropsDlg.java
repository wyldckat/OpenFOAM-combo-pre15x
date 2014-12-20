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

Description

\*---------------------------------------------------------------------------*/

package PatchTool.Shapes.EditPanels;

import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.Dimension;
import java.awt.event.*;
import java.awt.event.ActionEvent;
import java.awt.Graphics;
import java.awt.GraphicsConfiguration;
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import javax.media.j3d.*;
import javax.vecmath.*;

import PatchTool.PatchTool;

public class LinePropsDlg extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    public static final int OK = 0;
    public static final int CANCEL = 1;

    protected int state_;


    //--------------------------------------------------------------------------
    /** Creates new form LinePropsDlg */
    public LinePropsDlg(java.awt.Frame parent, Color3f colour)
    {
        super(parent, "Line Properties", true);   // non modal.

        state_ = -1;

        initComponents();

        linePropsPanel_.setColor3f(colour);

        show();
    }

    //--------------------------------------------------------------------------

    public int getState()
    {
        return state_;
    }

    //--------------------------------------------------------------------------

    public LinePropsPanel getLinePropsPanel()
    {
        return linePropsPanel_;
    }

    //--------------------------------------------------------------------------

    private void initComponents()
    {
        // This panel
        getContentPane().setLayout(new GridBagLayout());

        GridBagConstraints c = new GridBagConstraints();

        // Add insets around all to avoid clutter
        c.insets = new Insets(2, 2, 2, 2);

        // Anchor all same way
        c.anchor = GridBagConstraints.WEST;

        // Fill horizontal
        c.fill = GridBagConstraints.HORIZONTAL;


        // LineProps panel
        linePropsPanel_ = new LinePropsPanel();
        PatchTool.setConstraint(0, 0, 1, 1, 1.0, 1.0, c);
        getContentPane().add(linePropsPanel_, c);

        // Ok button
        okButton_ = new javax.swing.JButton("Ok");
        okButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    state_ = OK;
                    setVisible(false);
                    dispose();
                }
            }
        );
        PatchTool.setConstraint(0, 1, 1, 1, 1.0, 1.0, c);
        getContentPane().add(okButton_, c);

        // Cancel button
        cancelButton_ = new javax.swing.JButton("Cancel");
        cancelButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    state_ = CANCEL;
                    setVisible(false);
                    dispose();
                }
            }
        );
        PatchTool.setConstraint(1, 1, 1, 1, 1.0, 1.0, c);
        getContentPane().add(cancelButton_, c);

        pack();
    }


    //--------------------------------------------------------------------------
    //---- UI components
    //--------------------------------------------------------------------------

    private LinePropsPanel linePropsPanel_;
    private javax.swing.JButton okButton_;
    private javax.swing.JButton cancelButton_;

    //--------------------------------------------------------------------------
}


