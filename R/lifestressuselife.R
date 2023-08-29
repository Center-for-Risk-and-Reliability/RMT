# Use-Life Calculator
# Developed by Dr. Reuel Smith, 2021-2022

lifestress.uselife <-function(LSQest,MLEest,SUse){
  UselifeLSQ <- life(LSQest,SUse)
  UselifeMLE <- life(MLEest,SUse)
  return(list(LSQest,MLEest))
}
