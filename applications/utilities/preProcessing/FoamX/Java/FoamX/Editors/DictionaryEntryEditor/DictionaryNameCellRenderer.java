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

package FoamX.Editors.DictionaryEntryEditor;

import javax.swing.JTable;
import javax.swing.table.*;
import java.awt.Component;

import FoamX.App;

public class DictionaryNameCellRenderer
    extends DefaultTableCellRenderer
{
    //--------------------------------------------------------------------------

    protected java.awt.Font selectedFont_;
    protected java.awt.Color selectedColor_;

    //--------------------------------------------------------------------------
    /** DictionaryNameCellRenderer constructor. */
    public DictionaryNameCellRenderer()
    {
        selectedFont_= new java.awt.Font("Dialog", java.awt.Font.BOLD, 10);
        selectedColor_= new java.awt.Color(173, 255, 47);
    }

    //--------------------------------------------------------------------------

    public Component getTableCellRendererComponent
    (
        JTable  table,
        Object  value,
        boolean isSelected,
        boolean hasFocus,
        int     row,
        int     column
    )
    {
        Component comp = super.getTableCellRendererComponent
        (
            table,
            value,
            isSelected,
            hasFocus,
            row,
            column
        );

        if (isSelected)
        {
            //comp.setFont(selectedFont_);
            comp.setBackground(selectedColor_);
        }
        else
        {
            comp.setBackground(java.awt.Color.white);
        }


        return comp;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}

