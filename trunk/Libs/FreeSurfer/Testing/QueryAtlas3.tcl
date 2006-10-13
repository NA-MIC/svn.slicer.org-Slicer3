# assume that the model /projects/birn/freesurfer/data/bert/surf/lh.pial has been read in as the first model, fourth after the default ones

proc QueryAtlasInit { {filename ""} } {
  
  # find the data
  set ::QA(filename) ""
  if { $filename != "" } {
    set ::QA(filename) $filename
  } else {
    set candidates {
      /projects/birn/freesurfer/data/bert/surf/lh.pial
      i:/fBIRN-AHM2006/fbph2-000670986943/surf/lh.pial
      c:/data/fBIRN-AHM2006/fbph2-000648622547/surf/lh.pial
    }
    foreach c $candidates {
      if { [file exists $c] } {
        set ::QA(filename) $c
        set ::QA(directory) [file dirname [file dirname $::QA(filename)]]
        break
      }
    }
  }

  QueryAtlasAddModel
  QueryAtlasAddVolumes
  QueryAtlasAddAnnotations 
  QueryAtlasInitializePicker 
  QueryAtlasRenderView

  QueryAtlasUpdateCursor
}

#
# Add the model with the filename to the scene
#
proc QueryAtlasAddModel {} {

  # load the data
  set modelNode [vtkMRMLModelNode New]
  set modelStorageNode [vtkMRMLModelStorageNode New]
  set modelDisplayNode [vtkMRMLModelDisplayNode New]

  $modelStorageNode SetFileName $::QA(filename)
  if { [$modelStorageNode ReadData $modelNode] != 0 } {
    $modelNode SetName [file tail $::QA(filename)]

    $modelNode SetScene $::slicer3::MRMLScene
    $modelStorageNode SetScene $::slicer3::MRMLScene
    $modelDisplayNode SetScene $::slicer3::MRMLScene

    $::slicer3::MRMLScene AddNode $modelStorageNode
    $::slicer3::MRMLScene AddNode $modelDisplayNode

    $modelNode SetStorageNodeID [$modelStorageNode GetID]
    $modelNode SetAndObserveDisplayNodeID [$modelDisplayNode GetID]

    $::slicer3::MRMLScene AddNode $modelNode
    set ::QA(modelNodeID) [$modelNode GetID]
  }

  $modelNode Delete
  $modelStorageNode Delete
  $modelDisplayNode Delete
}

proc QueryAtlasAddVolumes {} {

  set selectionNode [$::slicer3::ApplicationLogic GetSelectionNode]
  set volumesLogic [$::slicer3::VolumesGUI GetLogic]
  set centered 1

  #
  # add the brain image
  #
  set fileName $::QA(directory)/mri/brain.mgz

  set volumeNode [$volumesLogic AddArchetypeVolume $fileName $centered]

  set ::QA(brain,volumeNodeID) [$volumeNode GetID]

  set volumeDisplayNode [$volumeNode GetDisplayNode]

  $volumeDisplayNode SetWindow 216
  $volumeDisplayNode SetLevel 108
  $volumeDisplayNode SetUpperThreshold 216
  $volumeDisplayNode SetLowerThreshold 30.99
  $volumeDisplayNode SetApplyThreshold 1

  #
  # add the function image
  #
  set fileName [file dirname $::QA(directory)]/sirp-hp65-stc-to7-gam.feat/stats/zstat8.nii

  set volumeNode [$volumesLogic AddArchetypeVolume $fileName $centered]

  set ::QA(functional,volumeNodeID) [$volumeNode GetID]
  set volumeDisplayNode [$volumeNode GetDisplayNode]


  #
  # add the segmentation image
  #
  set fileName $::QA(directory)/mri/aparc+aseg.mgz

  set volumeNode [$volumesLogic AddArchetypeVolume $fileName $centered]
  $volumeNode SetLabelMap 1

  set ::QA(label,volumeNodeID) [$volumeNode GetID]

  set volumeDisplayNode [$volumeNode GetDisplayNode]

  #
  # make brain be background and segmentation be label map
  #
  $selectionNode SetActiveVolumeID $::QA(brain,volumeNodeID)
  $selectionNode SetActiveLabelVolumeID $::QA(label,volumeNodeID)
  $::slicer3::ApplicationLogic PropagateVolumeSelection
  
  #
  # TODO: set up the real freesurfer aparc+aseg color map
  # - needs the labels to not go through window/level
  #
  EditorSetRandomLabelColormap 2500
}

#
# use the freesurfer annotation code to put 
# label scalars onto the model
#
proc QueryAtlasAddAnnotations {} {

  set fileName $::QA(filename)/../../label/lh.aparc.annot

  # get the model out of the scene
  set modelNode [$::slicer3::MRMLScene GetNodeByID $::QA(modelNodeID)]
  set displayNodeID [$modelNode GetDisplayNodeID]
  set displayNode [$::slicer3::MRMLScene GetNodeByID $displayNodeID]
  set viewer [$::slicer3::ApplicationGUI GetViewerWidget] 
  $viewer UpdateFromMRML
  set actor [$viewer GetActorByID [$modelNode GetID]]
  set mapper [$actor GetMapper]

  if [file exists $fileName] {

    set polydata [$modelNode GetPolyData]
    set scalaridx [[$polydata GetPointData] SetActiveScalars "labels"]

    if { $scalaridx == "-1" } {
        set scalars [vtkIntArray New]
        $scalars SetName "labels"
        [$polydata GetPointData] AddArray $scalars
        [$polydata GetPointData] SetActiveScalars "labels"
        $scalars Delete
    } 
    set scalaridx [[$polydata GetPointData] SetActiveScalars "labels"]
    set scalars [[$polydata GetPointData] GetArray $scalaridx]

    set lut [vtkLookupTable New]
    set fssar [vtkFSSurfaceAnnotationReader New]

    $fssar SetFileName $fileName
    $fssar SetOutput $scalars
    $fssar SetColorTableOutput $lut
    # try reading an internal colour table first
    $fssar UseExternalColorTableFileOff

    set retval [$fssar ReadFSAnnotation]
    if {$retval == 6} {
        error "ERROR: no internal colour table"
    }

    # set the look up table
    $mapper SetLookupTable $lut
    
    array unset _labels
    array set _labels [$fssar GetColorTableNames]
    array unset ::vtkFreeSurferReadersLabels_$::QA(modelNodeID)
    array set ::vtkFreeSurferReadersLabels_$::QA(modelNodeID) [array get _labels]
    set entries [lsort -integer [array names _labels]]

    # print them out
    set ::QA(labelMap) [array get _labels]

    # make the scalars visible
    $mapper SetScalarRange  [lindex $entries 0] [lindex $entries end]
    $mapper SetScalarVisibility 1

    $lut Delete
    $fssar Delete
    [$viewer GetMainViewer] Reset
  }
}

#
# convert a number to an RGBA 
# - A is always 255 (on transp)
# - number is incremented first so that 0 means background
#
proc QueryAtlasNumberToRGBA {number} {
  set number [expr $number + 1]
  set r [expr $number / (256 * 256)]
  set number [expr $number % (256 * 256)]
  set g [expr $number / 256]
  set b [expr $number % 256]

  return "$r $g $b 255"
}

#
# convert a RGBA to number
# - decrement by 1 to avoid ambiguity, since 0 is background
#
proc QueryAtlasRGBAToNumber {rgba} {
  foreach {r g b a} $rgba {}
  return [expr $r * (256*256) + $g * 256 + $b - 1] 
}


#
# set up a picking version of the polyData that can be used
# to render to the back buffer
#
proc QueryAtlasInitializePicker {} {

  #
  # add a prop picker to figure out if the mouse is actually over the model
  # and to identify the slice planes
  # and a pickRenderer and picker to find the actual world space pick point
  # under the mouse for a slice plane
  #
  set ::QA(propPicker) [vtkPropPicker New]
  set ::QA(cellPicker) [vtkCellPicker New]

  #
  # get the polydata for the model
  # - model node comes from the scene (retrieved by the ID)
  # - actor comes from the main Viewer
  # - mapper comes from the actor
  #
  set modelNode [$::slicer3::MRMLScene GetNodeByID $::QA(modelNodeID)]
  set ::QA(polyData) [vtkPolyData New]
  $::QA(polyData) DeepCopy [$modelNode GetPolyData]
  set ::QA(actor) [vtkActor New]
  set ::QA(mapper) [vtkPolyDataMapper New]
  $::QA(mapper) SetInput $::QA(polyData)
  $::QA(actor) SetMapper $::QA(mapper)

  #
  # instrument the polydata with cell number colors
  # - note: even though the array is named CellNumberColors here,
  #   vtk will (sometimes?) rename it to "Opaque Colors" as part of the first 
  #   render pass
  #

  $::QA(polyData) Update

  set cellData [$::QA(polyData) GetCellData]
  set cellNumberColors [$cellData GetArray "CellNumberColors"] 
  if { $cellNumberColors == "" } {
    set cellNumberColors [vtkUnsignedCharArray New]
    $cellNumberColors SetName "CellNumberColors"
    $cellData AddArray $cellNumberColors
    $cellData SetScalars $cellNumberColors
  }
  $cellData SetScalars $cellNumberColors

  set cellNumberColors [$cellData GetArray "CellNumberColors"] 
  $cellNumberColors Initialize
  $cellNumberColors SetNumberOfComponents 4

  set numberOfCells [$::QA(polyData) GetNumberOfCells]
  for {set i 0} {$i < $numberOfCells} {incr i} {
    eval $cellNumberColors InsertNextTuple4 [QueryAtlasNumberToRGBA $i]
  }

  set ::QA(cellData) $cellData
  set ::QA(numberOfCells) $numberOfCells

  set scalarNames {"CellNumberColors" "Opaque Colors"}
  foreach scalarName $scalarNames {
    if { [$::QA(cellData) GetScalars $scalarName] != "" } {
      $::QA(cellData) SetActiveScalars $scalarName
      break
    }
  }
  $::QA(mapper) SetScalarModeToUseCellData
  $::QA(mapper) SetScalarVisibility 1
  $::QA(mapper) SetScalarMaterialModeToAmbient
  $::QA(mapper) SetScalarRange 0 $::QA(numberOfCells)
  [$::QA(actor) GetProperty] SetAmbient 1.0
  [$::QA(actor) GetProperty] SetDiffuse 0.0

  #
  # add the mouse move callback
  # - create the classes that will be used every render and in the callback
  # - add the callback with the current render info
  #
  if { ![info exists ::QA(viewer)] } {
    set ::QA(viewer) [vtkImageViewer New]
    set ::QA(windowToImage) [vtkWindowToImageFilter New]
  }
  set renderWidget [[$::slicer3::ApplicationGUI GetViewerWidget] GetMainViewer]
  set renderer [$renderWidget GetRenderer]
  set interactor [$renderWidget GetRenderWindowInteractor] 
  set style [$interactor GetInteractorStyle] 

  $interactor AddObserver EnterEvent "QueryAtlasCursorVisibility on"
  $interactor AddObserver LeaveEvent "QueryAtlasCursorVisibility off"
  $interactor AddObserver MouseMoveEvent "QueryAtlasPickCallback"
  $interactor AddObserver RightButtonPressEvent "QueryAtlasMenuCreate start"
  $interactor AddObserver RightButtonReleaseEvent "QueryAtlasMenuCreate end"
  $style AddObserver StartInteractionEvent "QueryAtlasCursorVisibility off"
  $style AddObserver EndInteractionEvent "QueryAtlasCursorVisibility on"
  $style AddObserver EndInteractionEvent "QueryAtlasRenderView"

  $renderer AddActor $::QA(actor)
  $::QA(actor) SetVisibility 1
}


#
# re-render the picking model from the current camera location
#
proc QueryAtlasRenderView {} {

  #
  # get the renderer related instances
  #
  set renderWidget [[$::slicer3::ApplicationGUI GetViewerWidget] GetMainViewer]
  set renderWindow [$renderWidget GetRenderWindow]
  set renderer [$renderWidget GetRenderer]

  #
  # draw the image and get the pixels
  # - set the render parameters to draw with the cell labels
  # - draw in the back buffer
  # - pull out the pixels
  # - restore the draw state and render
  #
  $renderWindow SetSwapBuffers 0
  set renderState [QueryAtlasOverrideRenderState $renderer]
  $renderWidget Render


  $::QA(windowToImage) SetInput [$renderWidget GetRenderWindow]

  if { [$renderWindow GetSize] != [$::QA(viewer) GetSize] } {
    $::QA(windowToImage) Delete
    set ::QA(windowToImage) [vtkWindowToImageFilter New]
    $::QA(viewer) Delete
    set ::QA(viewer) [vtkImageViewer New]
    $::QA(windowToImage) SetInput [$renderWidget GetRenderWindow]
  }

  $::QA(viewer) SetColorWindow 255
  $::QA(viewer) SetColorLevel 127.5
  $::QA(windowToImage) SetInputBufferTypeToRGBA
  $::QA(windowToImage) ShouldRerenderOn
  $::QA(windowToImage) ReadFrontBufferOff
  $::QA(windowToImage) Modified
  $::QA(viewer) SetInput [$::QA(windowToImage) GetOutput]
  $::QA(viewer) Render

  $renderWindow SetSwapBuffers 1
  QueryAtlasRestoreRenderState $renderer $renderState
  $renderWidget Render

}

#####################################
# Override/Restore render state 
# - set up for rendering the cell picker and then
#   restore the state afterwards
#

proc QueryAtlasOverrideRenderState {renderer} {

  #
  # save the render state before overriding it with the 
  # parameters needed for cell rendering
  # - is just background color and visibility state of all actors
  #


  set actors [$renderer GetActors]
  set numberOfItems [$actors GetNumberOfItems]
  for {set i 0} {$i < $numberOfItems} {incr i} {
    set actor [$actors GetItemAsObject $i]
    set state($i,visibility) [$actor GetVisibility]
    $actor SetVisibility 0
  }

  set state(background) [$renderer GetBackground]
  $renderer SetBackground 0 0 0
  $::QA(actor) SetVisibility 1

  return [array get state]
}

proc QueryAtlasRestoreRenderState {renderer renderState} {

  array set state $renderState

  eval $renderer SetBackground $state(background)

  set actors [$renderer GetActors]
  set numberOfItems [$actors GetNumberOfItems]
  for {set i 0} {$i < $numberOfItems} {incr i} {
    set actor [$actors GetItemAsObject $i]
    $actor SetVisibility $state($i,visibility)
  }
  $::QA(actor) SetVisibility 0
}

### utility routine that should be provided by vtkCell
proc QueryAtlasPCoordsToWorld {cell pCoords} {
    
  if { [$cell GetClassName] != "vtkQuad" } {
    return "0 0 0"
  }

  foreach {r s t} $pCoords {}
  set rm [expr 1. - $r]
  set sm [expr 1. - $s]
  set sf0 [expr $rm * $sm]
  set sf1 [expr $r * $sm]
  set sf2 [expr $r * $s]
  set sf3 [expr $rm * $s]

  set points [$cell GetPoints]
  foreach {x0 x1 x2} "0 0 0" {}
  foreach p "0 1 2 3" {
    set point [$points GetPoint $p]
    foreach c "0 1 2" pp $point {
      set x$c [expr [set x$c] + $pp * [set sf$p]]
    }
  }
  return "$x0 $x1 $x2"
}


#
# query the cell number at the mouse location
#
proc QueryAtlasPickCallback {} {

  if { ![info exists ::QA(viewer)] } {
    return
  }

  if { [$::QA(cursor,actor) GetVisibility] == 0 } {
    # if the cursor isn't on, don't bother to calculate labels
    return
  }

  #
  # get access to the standard view parts
  #
  set viewer [$::slicer3::ApplicationGUI GetViewerWidget] 
  set renderWidget [$viewer GetMainViewer]
  set renderWindow [$renderWidget GetRenderWindow]
  set interactor [$renderWidget GetRenderWindowInteractor] 
  set renderer [$renderWidget GetRenderer]
  set actor [$viewer GetActorByID $::QA(modelNodeID)]

  # if the window size has changed, re-render
  if { [$renderWindow GetSize] != [$::QA(viewer) GetSize] } {
    QueryAtlasRenderView
  }

  # 
  # get the event location
  #
  eval $interactor UpdateSize [$renderer GetSize]
  set ::QA(lastWindowXY) [$interactor GetEventPosition]
  foreach {x y} $::QA(lastWindowXY) {}
  set ::QA(lastRootXY) [winfo pointerxy [$renderWidget GetWidgetName]]


  #
  # use the prop picker to see if we're over the model, or the slices
  # - set the 'hit' variable accordingly for later processing
  #
  set hit ""
  if { [$::QA(propPicker) PickProp $x $y $renderer] } {
    set prop [$::QA(propPicker) GetViewProp]
    if { $prop == $actor} {
      set hit "QueryActor"
    } else {
      set id [$viewer GetIDByActor $prop]
      set hit "Model"
    }
  }


  #
  # set the 'pointlabels' depending on the thing picked
  #

  switch $hit {

    "Model" {

      set node [$::slicer3::MRMLScene GetNodeByID $id]
      if { [$node GetDescription] != "" } {
        array set nodes [$node GetDescription]
        set nodes(sliceNode) [$::slicer3::MRMLScene GetNodeByID $nodes(SliceID)]
        set nodes(compositeNode) [$::slicer3::MRMLScene GetNodeByID $nodes(CompositeID)]

        set propCollection [$::QA(cellPicker) GetPickList]
        $propCollection RemoveAllItems
        $propCollection AddItem [$::QA(propPicker) GetViewProp]
        $::QA(cellPicker) PickFromListOn
        $::QA(cellPicker) Pick $x $y 0 $renderer
        set cellID [$::QA(cellPicker) GetCellId]
        set pCoords [$::QA(cellPicker) GetPCoords]
        if { $cellID != -1 } {
          set polyData [[$prop GetMapper] GetInput]
          set cell [$polyData GetCell $cellID]
          set rasPoint [QueryAtlasPCoordsToWorld $cell $pCoords]
          set pointLabels $rasPoint
        }
      }

    }

    "QueryActor" {

      #
      # get the color under the mouse from label image
      #
      $::QA(windowToImage) Update
      set color ""
      foreach c {0 1 2 3} {
        lappend color [[$::QA(windowToImage) GetOutput] GetScalarComponentAsFloat $x $y 0 $c]
      }

      #
      # convert the color to a cell index and get the cooresponding
      # label names from the vertices
      #
      set cellNumber [QueryAtlasRGBAToNumber $color]
      if { $cellNumber >= 0 && $cellNumber < [$::QA(polyData) GetNumberOfCells] } {
        set cell [$::QA(polyData) GetCell $cellNumber]

        set labels [[$::QA(polyData) GetPointData] GetScalars "labels"]

        array set labelMap $::QA(labelMap)
        set pointLabels ""
        set numberOfPoints [$cell GetNumberOfPoints]

        for {set p 0} {$p < $numberOfPoints} {incr p} {
          set index [$cell GetPointId $p]
          set pointLabel [$labels GetValue $index]
          if { [info exists labelMap($pointLabel)] } {
            set labelName $labelMap($pointLabel)
            if { [lsearch $pointLabels $labelName] == -1 } {
              lappend pointLabels $labelName
            }
          } else {
            lappend pointLabels "unknown"
          }
        }
        regsub -all " " $pointLabels "/" pointLabels
      } else {
        set pointLabels "background"
      }
    }
  }

  #
  # update the label
  #
  if { ![info exists ::QA(lastLabels)] } {
    set ::QA(lastLabels) ""
  }
  
  if { $pointLabels != $::QA(lastLabels) } {
    set ::QA(lastLabels) $pointLabels 
  }

  QueryAtlasUpdateCursor
}

proc QueryAtlasCursorVisibility { onoff } {

  if { $onoff == "on" } {
    $::QA(cursor,actor) SetVisibility 1
  } else {
    $::QA(cursor,actor) SetVisibility 0
  }
  set viewer [$::slicer3::ApplicationGUI GetViewerWidget]
  $viewer RequestRender

}

proc QueryAtlasUpdateCursor {} {

  set viewer [$::slicer3::ApplicationGUI GetViewerWidget]
  if { ![info exists ::QA(cursor,actor)] } {

    set ::QA(cursor,actor) [vtkTextActor New]
    set ::QA(cursor,mapper) [vtkTextMapper New]
    $::QA(cursor,actor) SetMapper $::QA(cursor,mapper)

    set renderWidget [$viewer GetMainViewer]
    set renderer [$renderWidget GetRenderer]
    $renderer AddActor2D $::QA(cursor,actor)

  }

  $::QA(cursor,actor) SetInput $::QA(lastLabels) 
  eval $::QA(cursor,actor) SetPosition $::QA(lastWindowXY) 
  $viewer RequestRender
}

proc QueryAtlasMenuCreate { state } {

  set renderWidget [[$::slicer3::ApplicationGUI GetViewerWidget] GetMainViewer]
  set interactor [$renderWidget GetRenderWindowInteractor] 
  set position [$interactor GetEventPosition]

  #
  # save the event position when the menu action started (when the right mouse
  # button was pressed) and only post the menu if the position is the same.  
  # If they aren't the same, do nothing since this was a dolly(zoom) action.
  #
  switch $state {
    "start" {
      set ::QA(menu,startPosition) $position
    }
    "end" {
      if { $::QA(menu,startPosition) == $position } {
        set parent [[$::slicer3::ApplicationGUI GetMainSlicerWindow] GetWidgetName]
        set qaMenu $parent.qaMenu
        catch "destroy $qaMenu"

        menu $qaMenu
        $qaMenu insert end command -label "Google..." -command "QueryAtlasQuery google"
        $qaMenu insert end command -label "Wikipedia..." -command "QueryAtlasQuery wikipedia"

        
        eval $qaMenu post $::QA(lastRootXY)
      }
    }
  }
}


proc QueryAtlasQuery { site } {

  switch $site {
    "google" {
      $::slicer3::Application OpenLink http://www.google.com/search?q=$::QA(lastLabels)
    }
    "wikipedia" {
      $::slicer3::Application OpenLink http://www.google.com/search?q=$::QA(lastLabels)+site:en.wikipedia.org
    }
  }
}

