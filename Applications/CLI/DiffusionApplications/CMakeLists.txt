PROJECT( DiffusionFiltering )

SET(filters
  jointLMMSE
  dwiNoiseFilter
  dwiUNLM
  #StochasticTractographyFilter
  SlicerFinslerTractography
  ResampleDTI
  )
FOREACH(filter ${filters})
  ADD_SUBDIRECTORY(${filter})
ENDFOREACH(filter)
