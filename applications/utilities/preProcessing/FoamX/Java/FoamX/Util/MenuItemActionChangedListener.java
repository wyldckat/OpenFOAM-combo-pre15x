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

package FoamX.Util;

import java.beans.*;
import javax.swing.*;

public class MenuItemActionChangedListener
    implements PropertyChangeListener
{
    private JMenuItem menuItem_;

    public MenuItemActionChangedListener(JMenuItem menuItem)
    {
        super();
        menuItem_ = menuItem;
    }

    public void propertyChange(PropertyChangeEvent e)
    {
        String propertyName = e.getPropertyName();
        if (e.getPropertyName().equals(Action.NAME))
        {
            // Change the label on the menu item if the action name has changed.
            String text = (String)e.getNewValue();
            menuItem_.setText(text);
        }
        else if (propertyName.equals("enabled"))
        {
            // Disable the menu item if the action has been disabled.
            Boolean enabledState = (Boolean)e.getNewValue();
            menuItem_.setEnabled(enabledState.booleanValue());
        }
    }
}




