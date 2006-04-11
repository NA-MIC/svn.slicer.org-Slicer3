# Load the KWWidgets package

package require kwwidgets

# Process some command-line arguments
# The --test option here is used to run this example as a non-interactive test
# for software quality purposes. Ignore this feature in your own application.

set option_test [expr [lsearch -exact $argv "--test"] == -1 ? 0 : 1]

# Create the application
# If --test was provided, ignore all registry settings, and exit silently
# Restore the settings that have been saved to the registry, like
# the geometry of the user interface so far.

set app [vtkKWApplication New]
$app SetName "MRML Scene Example"
if {$option_test} {
  $app SetRegistryLevel 0
  $app PromptBeforeExitOff
}
$app RestoreApplicationSettingsFromRegistry

# Set a help link. Can be a remote link (URL), or a local file

$app SetHelpDialogStartingPage "http://www.kwwidgets.org"

# Add a window
# Set 'SupportHelp' to automatically add a menu entry for the help link

set win [vtkKWWindow New]
$win SupportHelpOn
$app AddWindow $win
$win Create


#################################3
# create the mrml scene
set ::scene [vtkMRMLScene New]

set index [[$win GetEditMenu] AddCommand "&Undo" "" "$::scene Undo; mrmlUpdateUndoRedoButtons"]
[$win GetEditMenu] SetItemAccelerator $index "Ctrl+Z"
set index [[$win GetEditMenu] AddCommand "&Redo" "" "$::scene Redo; mrmlUpdateUndoRedoButtons"]
[$win GetEditMenu] SetItemAccelerator $index "Ctrl+Y"


# for use when cutting and pasting into the console:
# set scenefile c:/pieper/bwh/slicer3/latest/Slicer3/Applications/GUI/Testing/mrmlScene.xml

set scenefile [file dirname [info script]]/mrmlScene.xml
$::scene SetURL $scenefile
$::scene Connect


#
# add the needed pieces to display the volumes
# - a slice node to specify the plane
# - a SliceComposite node to identify the volumes
# - a SliceLogic to generate the image data
# - a SliceLayerLogic for the reslicing and display
#

# a SliceLayerLogic  - background and foreground
set ::slicebgl [vtkSlicerSliceLayerLogic New]
$::slicebgl SetMRMLScene $::scene
set ::slicefgl [vtkSlicerSliceLayerLogic New]
$::slicefgl SetMRMLScene $::scene

# a slice node to be controlled by the slicecontrol
set slicen [vtkMRMLSliceNode New]
$::scene AddNode $slicen
#$::slicen SetOrientationToCoronal
$::slicen SetOrientationToAxial
$::slicen SetDimensions 128 128 1
$::slicen SetFieldOfView 128 128 128
$::slicen UpdateMatrices

# a SliceComposite node for the SliceLogic
set slicecn [vtkMRMLSliceCompositeNode New]
$::slicecn SetBackgroundVolumeID vol3
$::slicecn SetForegroundVolumeID vol1
$::scene AddNode $slicecn

# a SliceLogic 
set ::slicel [vtkSlicerSliceLogic New]
$::slicel SetMRMLScene $::scene
$::slicel SetBackgroundLayer $::slicebgl
$::slicebgl UpdateTransforms
$::slicel SetForegroundLayer $::slicefgl
$::slicefgl UpdateTransforms
$::slicel SetSliceCompositeNode $::slicecn
$::slicel SetSliceNode $::slicen

##############
#
# key matrices:
#

set ::matrices ""

proc mrmlMatrices {} {

    if { $::matrices != "" } {
        return
    }

    set ::matrices [vtkKWTopLevel New]
    $::matrices SetApplication $::app
    $::matrices Create
    $::matrices SetTitle "Matrices"
    $::matrices Display
    $::matrices SetGeometry 240x1100+730+15

    $::matrices SetDeleteWindowProtocolCommand "" "set ::matrices {}"

    set ::XYToSliceFrame [vtkKWFrameWithLabel New]
    $::XYToSliceFrame SetParent $::matrices
    $::XYToSliceFrame Create
    $::XYToSliceFrame SetLabelText "XYToSlice"
    pack [$::XYToSliceFrame GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::XYToSlice [vtkKWMatrix4x4 New]
    $::XYToSlice SetParent [$::XYToSliceFrame GetFrame]
    $::XYToSlice Create
    $::XYToSlice SetMatrix4x4 [$::slicen GetXYToSlice]
    pack [$::XYToSlice GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::SliceToRASFrame [vtkKWFrameWithLabel New]
    $::SliceToRASFrame SetParent $::matrices
    $::SliceToRASFrame Create
    $::SliceToRASFrame SetLabelText "SliceToRAS"
    pack [$::SliceToRASFrame GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::SliceToRAS [vtkKWMatrix4x4 New]
    $::SliceToRAS SetParent [$::SliceToRASFrame GetFrame]
    $::SliceToRAS Create
    $::SliceToRAS SetMatrix4x4 [$::slicen GetSliceToRAS]
    pack [$::SliceToRAS GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::XYToRASFrame [vtkKWFrameWithLabel New]
    $::XYToRASFrame SetParent $::matrices
    $::XYToRASFrame Create
    $::XYToRASFrame SetLabelText "XYToRAS"
    pack [$::XYToRASFrame GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::XYToRAS [vtkKWMatrix4x4 New]
    $::XYToRAS SetParent [$::XYToRASFrame GetFrame]
    $::XYToRAS Create
    $::XYToRAS SetMatrix4x4 [$::slicen GetXYToRAS]
    pack [$::XYToRAS GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::RASToIJKMatrixbg [vtkMatrix4x4 New]

    set ::RASToIJKFrame [vtkKWFrameWithLabel New]
    $::RASToIJKFrame SetParent $::matrices
    $::RASToIJKFrame Create
    $::RASToIJKFrame SetLabelText "Bg RASToIJK"
    pack [$::RASToIJKFrame GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::RASToIJK [vtkKWMatrix4x4 New]
    $::RASToIJK SetParent [$::RASToIJKFrame GetFrame]
    $::RASToIJK Create
    $::RASToIJK SetMatrix4x4 $::RASToIJKMatrixbg
    pack [$::RASToIJK GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::XYToIJKFrame [vtkKWFrameWithLabel New]
    $::XYToIJKFrame SetParent $::matrices
    $::XYToIJKFrame Create
    $::XYToIJKFrame SetLabelText "Bg XYToIJK"
    pack [$::XYToIJKFrame GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::XYToIJK [vtkKWMatrix4x4 New]
    $::XYToIJK SetParent [$::XYToIJKFrame GetFrame]
    $::XYToIJK Create
    $::XYToIJK SetMatrix4x4 [[$::slicebgl GetXYToIJKTransform] GetMatrix]
    pack [$::XYToIJK GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2

    set ::matricesDismiss [vtkKWPushButton New]
    $::matricesDismiss SetParent $::matrices
    $::matricesDismiss Create
    $::matricesDismiss SetText "Dismiss"
    $::matricesDismiss SetCommand "" "$::matrices Delete; set ::matrices {}"
    pack [$::matricesDismiss GetWidgetName] -side top -anchor nw -expand true -fill both -padx 2 -pady 2
    
}


proc mrmlUpdateMatrices {} {

    if { $::matrices == "" } {
        return
    }
    
    $::XYToSlice UpdateWidget
    $::SliceToRAS UpdateWidget
    $::XYToRAS UpdateWidget
    [$::slicebgl GetVolumeNode] GetRASToIJKMatrix $::RASToIJKMatrixbg
    $::RASToIJK UpdateWidget
    $::XYToIJK UpdateWidget
    return

    set ::nodeIJKToRASbg [vtkMatrix4x4 New]
    [$::slicebgl GetVolumeNode] GetIJKToRASMatrix $::nodeIJKToRAS 
    puts "--------------IJKToRAS from the BG Volume Node"
    puts [$::nodeIJKToRAS Print]

    set ::nodeRASToIJKbg [vtkMatrix4x4 New]
    [$::slicebgl GetVolumeNode] GetRASToIJKMatrix $::nodeRASToIJK 
    puts "--------------RASToIJK from the BG Volume Node"
    puts [$::nodeRASToIJK Print]


    if {0} {
        [$::slicefgl GetVolumeNode] GetIJKToRASMatrix $nodeIJKToRAS 
        puts "--------------IJKToRAS from the FG Volume Node"
        puts [$nodeIJKToRAS Print]

        [$::slicefgl GetVolumeNode] GetRASToIJKMatrix $nodeRASToIJK 
        puts "--------------RASToIJK from the FG Volume Node"
        puts [$nodeRASToIJK Print]
    }

    puts "--------------XYToSlice from the Slice Node"
    puts [[$::slicen GetXYToSlice] Print]

    puts "--------------SliceRAS from the Slice Node"
    puts [[$::slicen GetSliceToRAS] Print]


    puts "--------------XYToRAS from the Slice Node"
    puts [[$::slicen GetXYToRAS] Print]

    puts "--------------XYToIJK from the BG layer logic"
    puts {[$::slicebgl GetXYToIJKTransform]}
    puts [[$::slicebgl GetXYToIJKTransform] Print]
}

##### #####


#
# set up the image viewer to render
#

set renderwidget [vtkKWRenderWidget New]
$renderwidget SetParent $win
$renderwidget Create
pack [$renderwidget GetWidgetName] -side top -fill both -expand y -padx 0 -pady 0

# Create an image viewer
# Use the render window and renderer of the renderwidget

set viewer [vtkImageViewer2 New]
$viewer SetColorWindow 200
$viewer SetColorLevel 100
$viewer SetRenderWindow [$renderwidget GetRenderWindow] 
$viewer SetRenderer [$renderwidget GetRenderer] 
$viewer SetInput [$::slicel GetImageData]
$viewer SetupInteractor [[$renderwidget GetRenderWindow] GetInteractor]

$renderwidget ResetCamera

## button to bring up matrices view

set ::matricesButton [vtkKWPushButton New]
$::matricesButton SetParent $::win
$::matricesButton Create
$::matricesButton SetText "Matrices"
$::matricesButton SetCommand "" "mrmlMatrices; mrmlUpdateMatrices"
pack [$::matricesButton GetWidgetName] -side top -anchor nw -expand true -fill x -padx 2 -pady 2

# SliceControl widget
set ::slicec [vtkSlicerSliceControlGUI New]
$::slicec SetParent $win
$::slicec Create
$::slicec SetSliceNode $slicen
$::slicec SetMRMLScene $::scene
pack [$::slicec GetWidgetName] -side top -anchor nw -expand false -fill x -padx 2 -pady 2
$::slicec AddObserver ModifiedEvent mrmlUpdateUndoRedoButtons
$::slicec AddObserver ModifiedEvent mrmlRender


# TODO: orientation doesn't work yet
[$::slicec GetOrientationMenu] SetEnabled 0

# VolumeSelect widget
set ::volsel [vtkSlicerVolumeSelectGUI New]
$::volsel SetParent $win
$::volsel Create
$::volsel SetMRMLScene $::scene
$::volsel SetLabelText "Volume Select: "
pack [$::volsel GetWidgetName] -side top -anchor nw -expand false -fill x -padx 2 -pady 2
$::volsel SetBalloonHelpString \
"select a volume from the current mrml scene."

# Undo button
set ::Undopb [vtkKWPushButton New]
$::Undopb SetParent $win
$::Undopb Create
$::Undopb SetText "Undo (0)"
$::Undopb SetEnabled 0
$::Undopb SetAnchorToWest
pack [$::Undopb GetWidgetName] -side top -anchor nw -expand false -fill x -padx 2 -pady 2
$::Undopb SetBalloonHelpString \
"Undo manipulations of mrml scene."
$::Undopb SetCommand "" "$::scene Undo; mrmlUpdateUndoRedoButtons"

# Redo button
set ::Redopb [vtkKWPushButton New]
$::Redopb SetParent $win
$::Redopb Create
$::Redopb SetText "Redo (0)"
$::Redopb SetEnabled 0
$::Redopb SetAnchorToWest
pack [$::Redopb GetWidgetName] -side top -anchor nw -expand false -fill x -padx 2 -pady 2
$::Redopb SetBalloonHelpString \
"Redo manipulations of mrml scene."
$::Redopb SetCommand "" "$::scene Redo; mrmlUpdateUndoRedoButtons"

# Add a button for each mrml node type

set nodes "vtkMRMLScalarVolumeNode vtkMRMLModelNode"
foreach n $nodes {

  set pb [vtkKWPushButton New]
  $pb SetParent $win
  $pb Create
  $pb SetText "$n"
  $pb SetAnchorToWest
  pack [$pb GetWidgetName] -side top -anchor nw -expand false -fill x -padx 2 -pady 2
  $pb SetBalloonHelpString \
    "Add a $n to the current mrml scene."
  $pb SetCommand "" "mrmlAddNode $n"

}

proc mrmlAddNode {nodetype} {

    set n [$nodetype New]
    $n SetName Node-[clock seconds]
    $::scene SaveStateForUndo
    $::scene AddNode $n
    mrmlUpdateUndoRedoButtons
}



# -----------------------------------------------------------------------

# Create a multi-column list

set mcl [vtkKWMultiColumnListWithScrollbars New]
$mcl SetParent $win
$mcl Create
$mcl SetBalloonHelpString \
"A simple multicolumn list. Columns can be resized moved and sorted.\
Double-click on some entries to edit them."
set mclwidget [$mcl GetWidget]
$mclwidget MovableColumnsOn
$mcl SetWidth 0
$mclwidget SetPotentialCellColorsChangedCommand \
    $mclwidget "ScheduleRefreshColorsOfAllCellsWithWindowCommand"
$mclwidget SetColumnSortedCommand \
    $mclwidget "ScheduleRefreshColorsOfAllCellsWithWindowCommand"

# Add the columns make some of them editable

set col_index [$mclwidget AddColumn "Node Type"] 

set col_index [$mclwidget AddColumn "Name"] 
$mclwidget SetColumnAlignmentToCenter $col_index
$mclwidget ColumnEditableOn $col_index

set col_index [$mclwidget AddColumn "ID"] 
$mclwidget SetColumnAlignmentToCenter $col_index


pack [$mcl GetWidgetName] -side top -anchor nw -expand n -fill both -padx 2 -pady 2

#
# set the command to execute when a cell is edited
#

$mclwidget SetCellUpdatedCommand "" "mrmlUpdateName $mclwidget"

proc mrmlUpdateName {mcl row col txt} {

    set id [$mcl GetCellText $row 2]
    set name [$mcl GetCellText $row 1]
    set node [$::scene GetNodeByID $id]
    $::scene SaveStateForUndo $node
    $node SetName $name
    $::scene Modified
    mrmlUpdateUndoRedoButtons
}

# 
# set an observer to refresh the list box when the scene changes
#

$::scene AddObserver ModifiedEvent "mrmlFillMCL $mclwidget"

proc mrmlFillMCL {mcl} {
    $mcl DeleteAllRows
    set nitems [$::scene GetNumberOfNodes]
    for {set i 0} {$i < $nitems} {incr i} {
        set n [$::scene GetNthNode $i]
        $mcl InsertCellText $i 0 [$n GetClassName]
        $mcl InsertCellText $i 1 [$n GetName]
        $mcl InsertCellText $i 2 [$n GetID]
    }
}

# initialize with the current scene
mrmlFillMCL $mclwidget

proc mrmlUpdateUndoRedoButtons {} {

    set nundo [$::scene GetNumberOfUndoLevels]
    $::Undopb SetText "Undo ($nundo)"
    if { $nundo > 0 } {
        $::Undopb SetEnabled 1
    } else {
        $::Undopb SetEnabled 0
    }

    set nredo [$::scene GetNumberOfRedoLevels]
    $::Redopb SetText "Redo ($nredo)"
    if { $nredo > 0 } {
        $::Redopb SetEnabled 1
    } else {
        $::Redopb SetEnabled 0
    }
}

proc mrmlRender {} {

    # update the slice layers
    $::slicebgl UpdateTransforms
    $::slicefgl UpdateTransforms

    # make the image actor match the output of the reslicing
    set dimx [expr [lindex [$::slicen GetDimensions] 0] -1]
    set dimy [expr [lindex [$::slicen GetDimensions] 1] -1]
    [$::viewer GetImageActor] SetDisplayExtent 0 $dimx 0 $dimy 0 0

    # update the widgets and render
    mrmlUpdateMatrices
    $::viewer Render
}

# initialize with the current state
mrmlUpdateUndoRedoButtons
mrmlRender


# Start the application
# If --test was provided, do not enter the event loop and run this example
# as a non-interactive test for software quality purposes.


set ret 0
$win Display
if {!$option_test} {
  $app Start
  set ret [$app GetExitStatus]
}
$win Close

# Deallocate and exit

$win Delete
$app Delete

