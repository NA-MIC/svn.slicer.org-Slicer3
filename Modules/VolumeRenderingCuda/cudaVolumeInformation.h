#ifndef __CUDAVOLUMEINFORMATION_H__
#define __CUDAVOLUMEINFORMATION_H__

//! A datastructure located on the cudacard that holds the information of the volume.
typedef struct {
    void* SourceData;
    int InputDataType;
    
    int VolumeSizeX, VolumeSizeY, VolumeSizeZ;


    //! The Color Transfer Function with a size of FunctionSize * 3 (rgb)
    float* ColorTransferFunction;
    //! The Alpha Transfer Function with a size of FunctionSize
    float* AlphaTransferFunction;
    //! The Size of the above Function
    unsigned int FunctionSize;

    //! The minimum and Maximum Values of the Volume
    float MinValueX, MaxValueX,
          MinValueY, MaxValueY,
          MinValueZ, MaxValueZ;

    //! The minimal Threshold of the Input Color Value     
    int MinThreshold;
    //! The Maximum Threshold of the Input Color Value
    int MaxThreshold;

    //! The Voxel Sizes.
    float VoxelSizeX, VoxelSizeY, VoxelSizeZ;

    //! The stepping accuracy to raster along the ray.
    float SteppingSize;

} cudaVolumeInformation;

#endif /* __CUDAVOLUMEINFORMATION_H__ */
