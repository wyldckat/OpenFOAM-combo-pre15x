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

Class
    InteractiveNode

Description
    Interface for 3D shapes that can be switched on/off and removed

\*---------------------------------------------------------------------------*/

package PatchTool.InteractiveNodes;

import javax.media.j3d.*;
import com.sun.j3d.utils.picking.PickResult;

public interface InteractiveNode
{
    public abstract String getName();

    public abstract BranchGroup getRoot();

    // Enable/disable
    public abstract void show();
    public abstract void hide();
    public abstract boolean isShowing();
    public abstract void remove();
    public abstract void edit();
}