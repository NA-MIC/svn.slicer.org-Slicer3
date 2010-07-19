#include <QApplication>

#include <qMRMLTreeWidget.h>
#include <qMRMLSceneModel.h>
#include <qMRMLSceneTransformModel.h>
#include <qMRMLSortFilterProxyModel.h>

#include <vtkMRMLScene.h>

#include <vtkTimerLog.h>

// STD includes
#include <stdlib.h>
#include <iostream>

int qMRMLTreeWidgetTest1( int argc, char * argv [] )
{
  QApplication app(argc, argv);
  std::cout << std::endl<< "***************************************" << std::endl;
  vtkTimerLog* timer = vtkTimerLog::New();
  vtkMRMLScene* scene = vtkMRMLScene::New();
  scene->SetURL("/home/julien/data/Slicer/SPL_PNL_Brain_Atlas2008/brain_atlas_2008.mrml");
  timer->StartTimer();
  scene->Import();
  timer->StopTimer();
  std::cout << std::endl << "Loaded: " << timer->GetElapsedTime() << std::endl;
  timer->StartTimer();
  scene->Delete();
  timer->StopTimer();
  std::cout << std::endl << "Deleted: " << timer->GetElapsedTime() << std::endl;

  std::cout << std::endl<< "***************************************" << std::endl;
  scene = vtkMRMLScene::New();
  qMRMLSceneModel   sceneModel;
  sceneModel.setMRMLScene(scene);
  scene->SetURL("/home/julien/data/Slicer/SPL_PNL_Brain_Atlas2008/brain_atlas_2008.mrml");
  timer->StartTimer();
  scene->Import();
  timer->StopTimer();
  std::cout << "qMRMLSceneModel Loaded: " << timer->GetElapsedTime() << std::endl;
  timer->StartTimer();
  scene->Delete();
  timer->StopTimer();
  std::cout << "qMRMLSceneModel Deleted: " << timer->GetElapsedTime() << std::endl;

  std::cout << std::endl<< "***************************************" << std::endl;
  scene = vtkMRMLScene::New();
  qMRMLSceneTransformModel   transformModel;
  transformModel.setMRMLScene(scene);
  scene->SetURL("/home/julien/data/Slicer/SPL_PNL_Brain_Atlas2008/brain_atlas_2008.mrml");
   timer->StartTimer();
  scene->Import();
  timer->StopTimer();
  std::cout << "qMRMLSceneTransformModel Loaded: " << timer->GetElapsedTime() << std::endl;
  timer->StartTimer();
  scene->Delete();
  timer->StopTimer();
  std::cout << "qMRMLSceneTransformModel Deleted: " << timer->GetElapsedTime() << std::endl;

  std::cout << std::endl<< "***************************************" << std::endl;
  scene = vtkMRMLScene::New();
  qMRMLSceneTransformModel   transformModel2;
  transformModel2.setMRMLScene(scene);
  qMRMLSortFilterProxyModel  sortModel;
  sortModel.setSourceModel(&transformModel2);
  scene->SetURL("/home/julien/data/Slicer/SPL_PNL_Brain_Atlas2008/brain_atlas_2008.mrml");
   timer->StartTimer();
  scene->Import();
  timer->StopTimer();
  std::cout << "qMRMLSceneTransformModel(+qMRMLSortFilterProxyModel) Loaded: " << timer->GetElapsedTime() << std::endl;
  timer->StartTimer();
  scene->Delete();
  timer->StopTimer();
  std::cout << "qMRMLSceneTransformModel(+qMRMLSortFilterProxyModel) Deleted: " << timer->GetElapsedTime() << std::endl;


  std::cout << std::endl<< "***************************************" << std::endl;
  scene = vtkMRMLScene::New();
  qMRMLTreeWidget   mrmlItem;
  mrmlItem.setMRMLScene(scene);
  scene->SetURL("/home/julien/data/Slicer/SPL_PNL_Brain_Atlas2008/brain_atlas_2008.mrml");
   timer->StartTimer();
  scene->Import();
  timer->StopTimer();
  std::cout << "qMRMLTreeWidget Loaded: " << timer->GetElapsedTime() << std::endl;
  timer->StartTimer();
  scene->Delete();
  timer->StopTimer();
  std::cout << "qMRMLTreeWidget Deleted: " << timer->GetElapsedTime() << std::endl;

  std::cout << std::endl<< "***************************************" << std::endl;
  scene = vtkMRMLScene::New();
  qMRMLTreeWidget   treeWidget;
  treeWidget.show();
  treeWidget.setMRMLScene(scene);
  scene->SetURL("/home/julien/data/Slicer/SPL_PNL_Brain_Atlas2008/brain_atlas_2008.mrml");
  timer->StartTimer();
  scene->Import();
  timer->StopTimer();
  std::cout << "qMRMLTreeWidget visible Loaded: " << timer->GetElapsedTime() << std::endl;
  timer->StartTimer();
  scene->Delete();
  timer->StopTimer();
  std::cout << "qMRMLTreeWidget visible Deleted: " << timer->GetElapsedTime() << std::endl;


  return EXIT_SUCCESS;
}
