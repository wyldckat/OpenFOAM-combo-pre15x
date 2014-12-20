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
package FoamX.Editors.DictionaryEntryEditor.EntryCache;

import java.text.*;
import java.awt.Frame;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import FoamXServer.CaseServer.*;
import FoamXServer.IDictionaryEntry;
import FoamXServer.IDictionaryEntryHolder;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.FoamXType;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.Util.FoamXAny;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;

public class LabelCache
    extends DictionaryEntryCache
//    implements DictionaryEntry
{

    //--------------------------------------------------------------------------
    /** LabelCache constructor. */
    public LabelCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** LabelCache constructor. */
    public LabelCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** LabelCache constructor. */
    public LabelCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------

    protected void initialiseEditor()
    {
        // Give dictionaryEntryCache opportunity to install combo editor
        super.initialiseEditor();

        if (editor_ == null)
        {
            // No custom editor - use string-based
            // DefaultCellEditor.
        }
    }

    //--------------------------------------------------------------------------
    //---- DictionaryEntry Interface Methods
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    /**
     *  Called by the DictionaryCellEditor's getCellEditorValue
     *  method to update the value of this item after editing.
     *  Will not be called by editors for compound types since e.g.
     *  dictionary elements cannot be added/removed in FoamX.
     */
    public boolean updateValue(Object value)
    {
        boolean bRet = false;
        String strValue = "";

        try
        {
            if (value != null)
            {
                // Incoming value object is a string.
                strValue = (String)value;

                // Check against type min and max if available
                FoamXAny minVal = typeDescriptor_.getMinValue();
                FoamXAny maxVal = typeDescriptor_.getMaxValue();
                if
                (
                    (
                        minVal.getAny().type.value()
                     == FoamXType._Type_Label
                    )
                 &&
                    (
                        maxVal.getAny().type.value()
                     == FoamXType._Type_Label
                    )
                )
                {
                    int intValue = Integer.parseInt(strValue);
                    int minValue = Integer.parseInt(minVal.toString());
                    int maxValue = Integer.parseInt(maxVal.toString());

                    if ((intValue < minValue) || (intValue > maxValue))
                    {
                        throw new FoamXException
                        (
                            "Label " + intValue + " out of range ["
                            + minValue + " , " + maxValue + "]"
                        );
                    }
                }

                // Update the Any object.
                value_.setValue(strValue);
                // Send the Any object to the dictionary entry object.
                setEntryValue();

                bRet = true;
            }
        }
        catch (NumberFormatException ex)
        {
            App.handleAllExceptions
            (
                new FoamXException
                (
                    "Number " + strValue + " not of correct type."
                )
            );
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return bRet;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




