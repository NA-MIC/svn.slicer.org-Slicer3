#ifndef __vtkSlicerFoundationIcons_h
#define __vtkSlicerFoundationIcons_h

#include "vtkKWObject.h"
#include "vtkKWResourceUtilities.h"
#include "vtkKWIcon.h"
#include "vtkSlicerIcons.h"
#include "./Resources/vtkSlicerFoundation_ImageData.h"

class VTK_SLICER_BASE_GUI_EXPORT vtkSlicerFoundationIcons : public vtkSlicerIcons
{
 public:
    static vtkSlicerFoundationIcons* New ( );
    vtkTypeRevisionMacro (vtkSlicerFoundationIcons, vtkSlicerIcons );
    void PrintSelf ( ostream& os, vtkIndent indent );

    vtkGetObjectMacro ( SlicerSelectAllIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerDeselectAllIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerTableIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerMoreOptionsIcon, vtkKWIcon);
    vtkGetObjectMacro ( SlicerGoIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing0Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing1Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing2Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing3Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing4Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing5Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing6Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoing7Icon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCameraIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerBlankIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCancelIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCancelDisabledIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCancelledIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCancelRequestedIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCleanUpIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerColorsIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerCycleIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerDecrementIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerDeleteIcon, vtkKWIcon);
    vtkGetObjectMacro ( SlicerDeleteDisabledIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerDoneIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerErrorIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGlyphIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerGoToFirstIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerGoToLastIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerIncrementIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerInformationIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerLoadIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerDownloadIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerUploadIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFoundOnDiskIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerNotFoundOnDiskIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerPauseIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerPingPongIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerBackwardIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerForwardIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPreparingIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerRecordIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerSaveIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPlayerStopRecordingIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerTimedOutIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerTinyHelpIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerWaitIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerMagnifyIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerMinifyIcon, vtkKWIcon );    
    vtkGetObjectMacro ( SlicerNextIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerPreviousIcon, vtkKWIcon );    
    vtkGetObjectMacro ( SlicerGoToEndIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerGoToStartIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerUndoIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerRedoIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerUnlinkIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerLinkIcon, vtkKWIcon );    
    vtkGetObjectMacro ( SlicerCheckedVisibleIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerVisibleNoFrameIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerInvisibleNoFrameIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerVisibleIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerInvisibleIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerRefreshIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerVolumeIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerCenterOnFiducialIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsAddNewIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsDeleteAllIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsDeleteAllInListIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsDeleteLastClickedIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsSelectAllIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsSelectNoneIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsSelectAllInListIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsSelectNoneInListIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsUpIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsDownIcon, vtkKWIcon );    
    vtkGetObjectMacro ( SlicerFiducialsDeleteListIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsLockListIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsUnlockListIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsHideListIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsExposeListIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerFiducialsHideExposeAllListsIcon, vtkKWIcon );    
    vtkGetObjectMacro (SlicerFiducialsHideAllListsIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerFiducialsExposeAllListsIcon, vtkKWIcon );    
    vtkGetObjectMacro (SlicerLockIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerUnlockIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerCompositeIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerLockOrUnlockIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerVisibleOrInvisibleIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerImportSceneIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadSceneIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadDicomVolumeIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadVolumeIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadDirectoryIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadFiducialsIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadModelIcon, vtkKWIcon);
    vtkGetObjectMacro (SlicerLoadTransformIcon, vtkKWIcon);
    vtkGetObjectMacro ( SlicerLoadColorLUTIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerLoadFiberBundleIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerLoadScalarOverlayIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerCloseIcon, vtkKWIcon );
    vtkGetObjectMacro ( SlicerExtensionsIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerWWWIcon, vtkKWIcon );
    vtkGetObjectMacro (SlicerRotateToPixelSpaceIcon, vtkKWIcon );

    void AssignImageDataToIcons ( );
    
 protected:
    vtkSlicerFoundationIcons ( );
    virtual ~vtkSlicerFoundationIcons ( );
    vtkKWIcon *SlicerLockOrUnlockIcon;
    vtkKWIcon *SlicerVisibleOrInvisibleIcon;
    vtkKWIcon *SlicerSelectAllIcon;
    vtkKWIcon *SlicerDeselectAllIcon;
    vtkKWIcon *SlicerTableIcon;
    vtkKWIcon *SlicerGoIcon;
    vtkKWIcon *SlicerGoing0Icon;
    vtkKWIcon *SlicerGoing1Icon;
    vtkKWIcon *SlicerGoing2Icon;
    vtkKWIcon *SlicerGoing3Icon;
    vtkKWIcon *SlicerGoing4Icon;
    vtkKWIcon *SlicerGoing5Icon;
    vtkKWIcon *SlicerGoing6Icon;
    vtkKWIcon *SlicerGoing7Icon;
    vtkKWIcon *SlicerCameraIcon;
    vtkKWIcon *SlicerBlankIcon;
    vtkKWIcon *SlicerCancelIcon;
    vtkKWIcon *SlicerCancelDisabledIcon;
    vtkKWIcon *SlicerCancelledIcon;
    vtkKWIcon *SlicerCancelRequestedIcon;
    vtkKWIcon *SlicerCleanUpIcon;
    vtkKWIcon *SlicerColorsIcon;
    vtkKWIcon *SlicerPlayerCycleIcon;
    vtkKWIcon *SlicerDecrementIcon;
    vtkKWIcon *SlicerDeleteIcon;
    vtkKWIcon *SlicerDeleteDisabledIcon;
    vtkKWIcon *SlicerDoneIcon;
    vtkKWIcon *SlicerErrorIcon;
    vtkKWIcon *SlicerGlyphIcon;
    vtkKWIcon *SlicerPlayerGoToFirstIcon;
    vtkKWIcon *SlicerPlayerGoToLastIcon;
    vtkKWIcon *SlicerIncrementIcon;
    vtkKWIcon *SlicerInformationIcon;
    vtkKWIcon *SlicerLoadIcon;
    vtkKWIcon *SlicerUploadIcon;
    vtkKWIcon *SlicerDownloadIcon;
    vtkKWIcon *SlicerFoundOnDiskIcon;
    vtkKWIcon *SlicerNotFoundOnDiskIcon;
    vtkKWIcon *SlicerPlayerPauseIcon;
    vtkKWIcon *SlicerPlayerPingPongIcon;
    vtkKWIcon *SlicerPlayerBackwardIcon;
    vtkKWIcon *SlicerPlayerForwardIcon;
    vtkKWIcon *SlicerPreparingIcon;
    vtkKWIcon *SlicerPlayerRecordIcon;
    vtkKWIcon *SlicerSaveIcon;
    vtkKWIcon *SlicerPlayerStopRecordingIcon;
    vtkKWIcon *SlicerTimedOutIcon;
    vtkKWIcon *SlicerTinyHelpIcon;
    vtkKWIcon *SlicerWaitIcon;
    vtkKWIcon *SlicerMagnifyIcon;
    vtkKWIcon *SlicerMinifyIcon;
    vtkKWIcon *SlicerNextIcon;
    vtkKWIcon *SlicerPreviousIcon;
    vtkKWIcon *SlicerGoToEndIcon;
    vtkKWIcon *SlicerGoToStartIcon;
    vtkKWIcon *SlicerUndoIcon;
    vtkKWIcon *SlicerRedoIcon;
    vtkKWIcon *SlicerUnlinkIcon;
    vtkKWIcon *SlicerLinkIcon;    
    vtkKWIcon *SlicerCheckedVisibleIcon;
    vtkKWIcon *SlicerVisibleNoFrameIcon;
    vtkKWIcon *SlicerInvisibleNoFrameIcon;
    vtkKWIcon *SlicerVisibleIcon;
    vtkKWIcon *SlicerInvisibleIcon;
    vtkKWIcon *SlicerRefreshIcon;
    vtkKWIcon *SlicerVolumeIcon;
    vtkKWIcon *SlicerMoreOptionsIcon;
    vtkKWIcon *SlicerCenterOnFiducialIcon;
    vtkKWIcon *SlicerFiducialsAddNewIcon;
    vtkKWIcon *SlicerFiducialsDeleteAllIcon;
    vtkKWIcon *SlicerFiducialsDeleteAllInListIcon;
    vtkKWIcon *SlicerFiducialsDeleteLastClickedIcon;
    vtkKWIcon *SlicerFiducialsSelectAllIcon;
    vtkKWIcon *SlicerFiducialsSelectNoneIcon;
    vtkKWIcon *SlicerFiducialsSelectAllInListIcon;
    vtkKWIcon *SlicerFiducialsSelectNoneInListIcon;
    vtkKWIcon *SlicerFiducialsUpIcon;
    vtkKWIcon *SlicerFiducialsDownIcon;
    vtkKWIcon *SlicerFiducialsDeleteListIcon;
    vtkKWIcon *SlicerFiducialsLockListIcon;
    vtkKWIcon *SlicerFiducialsUnlockListIcon;
    vtkKWIcon *SlicerFiducialsHideListIcon;
    vtkKWIcon *SlicerFiducialsExposeListIcon;
    vtkKWIcon *SlicerFiducialsHideExposeAllListsIcon;    
    vtkKWIcon *SlicerFiducialsHideAllListsIcon;
    vtkKWIcon *SlicerFiducialsExposeAllListsIcon;
    vtkKWIcon *SlicerCompositeIcon;
    vtkKWIcon *SlicerLockIcon;
    vtkKWIcon *SlicerUnlockIcon;
    vtkKWIcon *SlicerImportSceneIcon;
    vtkKWIcon *SlicerLoadSceneIcon;
    vtkKWIcon *SlicerLoadDicomVolumeIcon;
    vtkKWIcon *SlicerLoadVolumeIcon;
    vtkKWIcon *SlicerLoadDirectoryIcon;
    vtkKWIcon *SlicerLoadFiducialsIcon;
    vtkKWIcon *SlicerLoadModelIcon;
    vtkKWIcon *SlicerLoadTransformIcon;
    vtkKWIcon *SlicerLoadColorLUTIcon;
    vtkKWIcon *SlicerLoadFiberBundleIcon;
    vtkKWIcon *SlicerLoadScalarOverlayIcon;
    vtkKWIcon *SlicerCloseIcon;
    vtkKWIcon *SlicerExtensionsIcon;
    vtkKWIcon *SlicerWWWIcon;
    vtkKWIcon *SlicerRotateToPixelSpaceIcon;
    
 private:
    vtkSlicerFoundationIcons ( const vtkSlicerFoundationIcons&); /// Not implemented
    void operator = (const vtkSlicerFoundationIcons& ); /// not implemented.
    
};
#endif
