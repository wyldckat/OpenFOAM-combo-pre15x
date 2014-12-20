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
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundEditor;
import FoamX.CaseManagement.CaseChooserDlg;
import FoamX.CaseManagement.CaseBrowserPanel;

public class RootAndCaseCache
    extends DictionaryEntryCache
//    implements DictionaryEntry
{

    //--------------------------------------------------------------------------
    /** RootAndCaseCache constructor. */
    public RootAndCaseCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** RootAndCaseCache constructor. */
    public RootAndCaseCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** RootAndCaseCache constructor. */
    public RootAndCaseCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }


    /**
     * Display as <case>,<root>
     */
    public String toString()
    {
        try
        {
            displayValue_ = value_.toString();

            if (displayValue_.length() != 0)
            {
                int lastSlash = displayValue_.lastIndexOf('/');

                if (lastSlash == -1)
                {
                    throw new FoamXException
                    (
                        "RootAndCase " + displayValue_ + " cannot be split "
                        + " into root directory and case name"
                    );

                }
                displayValue_ =
                    displayValue_.substring(lastSlash+1)
                  + ','
                  + displayValue_.substring(0, lastSlash);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return displayValue_;
    }

    //--------------------------------------------------------------------------
    /**
     * For argument strings rootAndCase becomes space separated
     */
    public void toStringRaw(Vector argListVector)
    {
        try
        {
            String rootAndCase = value_.toString();

            if (rootAndCase.length() != 0)
            {
                int lastSlash = rootAndCase.lastIndexOf('/');

                if (lastSlash == -1)
                {
                    throw new FoamXException
                    (
                        "RootAndCase " + rootAndCase + " cannot be split "
                        + " into root directory and case name"
                    );

                }
                argListVector.addElement(rootAndCase.substring(0, lastSlash));
                argListVector.addElement(rootAndCase.substring(lastSlash+1));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }
    //--------------------------------------------------------------------------

    protected void initialiseEditor()
    {
        // Use compound editor.
        DictionaryEntryCompoundEditor compEdit =
            new DictionaryEntryCompoundEditor();

        // Subscribe to the edit button's action event so
        // that we can invoke the appropriate edit action for
        // the compound types.
        compEdit.addActionListener
        (
            new java.awt.event.ActionListener()
            {
                public void actionPerformed
                (
                    java.awt.event.ActionEvent evt
                )
                {
                    editButtonActionPerformed(evt);
                }
            }
        );
        editor_ = compEdit;
    }

    //--------------------------------------------------------------------------

    protected void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // The edit button has been pressed. The user wants to edit this entry.
        try
        {
            // Pop up case chooser panel
            CaseChooserDlg caseChooser =
                new CaseChooserDlg
                (
                    App.getRootFrame(),
                    CaseBrowserPanel.SELECT_CASE_MODE
                );
            caseChooser.show();

            if 
            (
                (caseChooser.getCaseRoot() != null)
             && (caseChooser.getCaseName() != null)
             && (caseChooser.getCaseBrowser() != null)
            )
            {
                updateValue
                (
                    caseChooser.getCaseRoot()
                  + '/'
                  + caseChooser.getCaseName()
                );

                // Signal that editing has stopped.
                editor_.stopCellEditing();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
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
        String strValue;

        try
        {
            if (value != null)
            {
                // Incoming value object is a string.
                strValue = (String)value;

                // Update the Any object.
                value_.setValue(strValue);
                // Send the Any object to the dictionary entry object.
                setEntryValue();

                bRet = true;
            }
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




