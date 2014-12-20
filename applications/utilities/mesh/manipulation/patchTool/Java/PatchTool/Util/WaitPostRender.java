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

Description

\*---------------------------------------------------------------------------*/

package PatchTool.Util;

import java.awt.GraphicsConfiguration;
import javax.media.j3d.*;
import javax.vecmath.*;
import java.awt.event.*;
import javax.swing.event.*;

import PatchTool.Events.RenderListener;

public class WaitPostRender implements RenderListener
{
    protected MyCanvas3D c_;

    //--------------------------------------------------------------------------

    public WaitPostRender(MyCanvas3D c)
    {
        c_ = c;

	c_.addRenderListener(this);
    }

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //---- RenderListener Implementation
    //--------------------------------------------------------------------------

    public void postSwap(java.util.EventObject evt)
    {
        c_.removeRenderListener(this);
        System.out.println("WaitPostRender : Done rendering");
    }

    //--------------------------------------------------------------------------

}
