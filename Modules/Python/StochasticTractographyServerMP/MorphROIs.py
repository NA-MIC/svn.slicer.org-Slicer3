import numpy as np
import scipy as si
import scipy.ndimage.morphology as mi  
import imp, glob, os, sys


def morph(nLayers, isDilate):
  struct = mi.generate_binary_structure(3,1).astype('uint16')

  print 'Current directory : ', os.getcwd()
  rois = glob.glob('*.roi')
  print 'Files : ', rois


  tmpD = 'morphed'
  isDir = os.access('morphed', os.F_OK)
  if not isDir:
     os.mkdir('morphed')  


  for i in range(len(rois)):
    dims = np.fromfile(rois[i].split('.')[0] + '.dims', 'uint16')
    data = np.fromfile(rois[i], 'uint16').reshape(dims[0], dims[1], dims[2])
    if isDilate:
      data = mi.binary_dilation(data, struct, nLayers).astype('uint16')
    else:
      data = mi.binary_erosion(data, struct, nLayers).astype('uint16')
    #data.tofile(rois[i].split('.')[0]+'_dilated' + '.roi')
    data.tofile(tmpD + '/' + rois[i].split('.')[0] + '_' + str(nLayers) + '_dilated'  + '.roi')

if __name__ == '__main__':
   if len(sys.argv)==3:
     nLayers = int(sys.argv[1])
     isDilate = bool(int(sys.argv[2]))
     print 'Number of layers : ', nLayers
     print 'Dilation : ', isDilate
   else:
     nLayers = 1
     isDilate = True

   morph(nLayers, isDilate)
