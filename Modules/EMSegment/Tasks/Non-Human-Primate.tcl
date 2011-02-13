package require Itcl

#########################################################
#
if {0} { ;# comment

    This is function is executed by EMSegmenter

    # TODO :

}
#
#########################################################

#
# namespace procs
#

#
# Remember to source first GenericTask.tcl as it has all the variables/basic structure defined
#
namespace eval EMSegmenterPreProcessingTcl {

    #
    # Variables Specific to this Preprocessing
    #
    variable TextLabelSize 1
    variable CheckButtonSize 2
    variable VolumeMenuButtonSize 0
    variable TextEntrySize 0

    # Check Button
    variable atlasAlignedFlagID 0
    variable inhomogeneityCorrectionFlagID 1

    # Text Entry
    # not defined for this task

    #
    # OVERWRITE DEFAULT
    #

    # -------------------------------------
    # Define GUI
    # return 1 when error occurs
    # -------------------------------------
    proc ShowUserInterface { } {
        variable preGUI
        variable atlasAlignedFlagID
        variable iccMaskSelectID
        variable inhomogeneityCorrectionFlagID
        variable LOGIC

        # Always has to be done initially so that variables are correctly defined
        if { [InitVariables] } {
            puts stderr "ERROR: MRI-HumanBrain: ShowUserInterface: Not all variables are correctly defined!"
            return 1
        }
        $LOGIC PrintText  "TCLMRI: Preprocessing MRI Human Brain - ShowUserInterface"

        $preGUI DefineTextLabel "This task only applies to non-skull stripped scans! \n\nShould the EMSegmenter " 0
        $preGUI DefineCheckButton "- register the atlas to the input scan ?" 0 $atlasAlignedFlagID
        # $preGUI DefineCheckButton "Are the input scans skull stripped ?" 0 $skullStrippedFlagID
        # $preGUI DefineVolumeMenuButton "Define ICC mask of the atlas ?" 0 $iccMaskSelectID
        $preGUI DefineCheckButton "- perform image inhomogeneity correction on input scan ?" 0 $inhomogeneityCorrectionFlagID

        # Define this at the end of the function so that values are set by corresponding MRML node
        $preGUI SetButtonsFromMRML
    }

    # -------------------------------------
    # Define Preprocessing Pipeline
    # return 1 when error occurs
    # -------------------------------------
    proc Run { } {
        variable preGUI
        variable workingDN
        variable alignedTargetNode
        variable inputAtlasNode
        variable mrmlManager
        variable LOGIC

        variable atlasAlignedFlagID
        variable iccMaskSelectID
        variable inhomogeneityCorrectionFlagID

        $LOGIC PrintText "TCLMRI: =========================================="
        $LOGIC PrintText "TCLMRI: == Preprocress Data"
        $LOGIC PrintText "TCLMRI: =========================================="
        # ---------------------------------------
        # Step 1 : Initialize/Check Input
        if {[InitPreProcessing]} { 
            return 1
        }

        # ----------------------------------------------------------------------------
        # We have to create this function so that we can run it in command line mode
        #
        set atlasAlignedFlag [ GetCheckButtonValueFromMRML $atlasAlignedFlagID ]
        set skullStrippedFlag 0
        set iccMaskVTKID 0
        # [GetVolumeMenuButtonValueFromMRML $iccMaskSelectID]
        set inhomogeneityCorrectionFlag [GetCheckButtonValueFromMRML $inhomogeneityCorrectionFlagID]

        $LOGIC PrintText "TCLMRI: ==> Preprocessing Setting: $atlasAlignedFlag $inhomogeneityCorrectionFlag"

        if { ($atlasAlignedFlag == 0) && ($skullStrippedFlag == 1) } {
            PrintError "Run: We currently cannot align the atlas to skull stripped image"
            return 1
        }

        if { $iccMaskVTKID } {
            set inputAtlasICCMaskNode [$mrmlManager GetVolumeNode $iccMaskVTKID]
            if { $inputAtlasICCMaskNode == "" } {
                PrintError "Run: inputAtlasICCMaskNode is not defined"
                return 1
            }
        } else {
            set inputAtlasICCMaskNode ""
        }

        # -------------------------------------
        # Step 2: Generate ICC Mask Of input images
        if { $inputAtlasICCMaskNode != "" && 0} {
            set inputAtlasVolumeNode [$inputAtlas GetNthVolumeNode 0]
            set alignedTargetVolumeNode [$alignedTargetNode GetNthVolumeNode 0]

            set subjectICCMaskNode [GenerateICCMask $inputAtlasVolumeNode $inputAtlasICCMaskNode $alignedTargetVolumeNode]

            if { $subjectICCMaskNode == "" } {
                PrintError "Run: Generating ICC mask for Input failed!"
                return 1
            }
        } else {
            #  $LOGIC PrintText "TCLMRI: Skipping ICC Mask generation! - Not yet implemented"
            set subjectICCMaskNode ""
        }

        # -------------------------------------
        # Step 4: Perform Intensity Correction
        if { $inhomogeneityCorrectionFlag == 1 } {

            set subjectIntensityCorrectedNodeList [PerformIntensityCorrection $subjectICCMaskNode]
            if { $subjectIntensityCorrectedNodeList == "" } {
                PrintError "Run: Intensity Correction failed !"
                return 1
            }
            if { [UpdateVolumeCollectionNode "$alignedTargetNode" "$subjectIntensityCorrectedNodeList"] } {
                return 1
            }
        } else {
             $LOGIC PrintText "TCLMRI: Skipping intensity correction"
        }

        # write results over to alignedTargetNode

        # -------------------------------------
        # Step 5: Atlas Alignment - you will also have to include the masks
        # Defines $workingDN GetAlignedAtlasNode
        if { [RegisterAtlas $atlasAlignedFlag] } {
            PrintError "Run: Atlas alignment failed !"
            return 1
        }


        # -------------------------------------
        # Step 6: Perform autosampling to define intensity distribution
        if { [ComputeIntensityDistributions] } {
            PrintError "Run: Could not automatically compute intensity distribution !"
            return 1
        }

        # -------------------------------------
        # Step 7: Check validity of Distributions 
        set failedIDList [CheckAndCorrectTreeCovarianceMatrix]
        if { $failedIDList != "" } {
            set MSG "Log Covariance matrices for the following classes seemed incorrect:\n "
            foreach ID $failedIDList {
                set MSG "${MSG}[$mrmlManager GetTreeNodeName $ID]\n"
            }
            set MSG "${MSG}This can cause failure of the automatic segmentation. To address the issue, please visit the web site listed under Help"

            if { 0 } {
                # TODO
                $preGUI PopUpWarningWindow "$MSG"
            } else {
                $LOGIC PrintText "TCLMRI: WARNING: $MSG"
            }
        }

        return 0
    }

    #
    # TASK SPECIFIC FUNCTIONS
    #

}


namespace eval EMSegmenterSimpleTcl {
    # 0 = Do not create a check list for the simple user interface
    # simply remove
    # 1 = Create one - then also define ShowCheckList and
    #     ValidateCheckList where results of checklist are transfered to Preprocessing

    proc CreateCheckList { } {
        return 1
    }

    proc ShowCheckList { } {
        variable inputChannelGUI
        # Always has to be done initially so that variables are correctly defined
        if { [InitVariables] } {
            PrintError "ShowCheckList: Not all variables are correctly defined!"
            return 1
        }

        $inputChannelGUI DefineTextLabel "Please insure that input scans are not skull stripped" 0
        $inputChannelGUI DefineCheckButton "Perform image inhomogeneity correction on input scans ?" 0 $EMSegmenterPreProcessingTcl::inhomogeneityCorrectionFlagID

        # Define this at the end of the function so that values are set by corresponding MRML node
        $inputChannelGUI SetButtonsFromMRML
        return 0

    }

    proc ValidateCheckList { } {
        return 0
    }

    proc PrintError { TEXT } {
         puts stderr "TCLMRI: ERROR:EMSegmenterSimpleTcl::${TEXT}"
    }
}

