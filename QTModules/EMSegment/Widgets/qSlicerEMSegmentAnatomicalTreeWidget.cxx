// Qt includes
#include <QDebug>
#include <QTreeView>
#include <QStandardItemModel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QToolButton>
#include <QCheckBox>

// CTK includes
#include <ctkLogger.h>

// EMSegment includes
#include "qSlicerEMSegmentAnatomicalTreeWidget.h"

// EMSegment/MRML includes
#include <vtkEMSegmentMRMLManager.h>
#include <vtkMRMLEMSTreeNode.h>

// MRML includes
#include <vtkMRMLVolumeNode.h>

//--------------------------------------------------------------------------
static ctkLogger logger(
    "org.slicer.qtmodules.emsegment.widgets.qSlicerEMSegmentAnatomicalTreeWidget");
//--------------------------------------------------------------------------

//-----------------------------------------------------------------------------
class qSlicerEMSegmentAnatomicalTreeWidgetPrivate :
    public ctkPrivate<qSlicerEMSegmentAnatomicalTreeWidget>
{
public:
  qSlicerEMSegmentAnatomicalTreeWidgetPrivate();

  enum
    {
    TreeNodeIDRole = Qt::UserRole + 1,
    TreeItemTypeRole
    };

  enum TreeItemType
    {
    StructureNameItemType = 0,
    LabelItemType,
    MRMLIDItemType,
    ClassWeightItemType,
    AtlasWeightItemType,
    AlphaItemType
    };

  QTreeView *              TreeView;
  QStandardItemModel *     TreeModel;
  bool                     StructureNameEditable;
  bool                     LabelColumnVisible;
  bool                     ClassWeightColumnVisible;
  bool                     AtlasWeightColumnVisible;
  bool                     AlphaColumnVisible;

  QCheckBox *              DisplayMRMLIDsCheckBox;
  QToolButton *            CollapseAllButton;
  QToolButton *            ExpandAllButton;
};

//-----------------------------------------------------------------------------
// qSlicerEMSegmentAnatomicalTreeWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerEMSegmentAnatomicalTreeWidgetPrivate::qSlicerEMSegmentAnatomicalTreeWidgetPrivate()
{
  this->TreeView = 0;
  this->TreeModel = new QStandardItemModel;

  this->StructureNameEditable = false;
  this->LabelColumnVisible = false;
  this->ClassWeightColumnVisible = false;
  this->AtlasWeightColumnVisible = false;
  this->AlphaColumnVisible = false;

  QStringList headerNames;
  headerNames << "Structure" << "Id" << "Label"
      << "Class Weight" << "Atlas Weight" << "Alpha";
  this->TreeModel->setHorizontalHeaderLabels(headerNames);
}

//-----------------------------------------------------------------------------
// qSlicerEMSegmentAnatomicalTreeWidget methods

//-----------------------------------------------------------------------------
qSlicerEMSegmentAnatomicalTreeWidget::qSlicerEMSegmentAnatomicalTreeWidget(QWidget *newParent):
Superclass(newParent)
{
  logger.setDebug();
  CTK_INIT_PRIVATE(qSlicerEMSegmentAnatomicalTreeWidget);
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);

  // Layout (TreeView and (control buttons)) and DisplayMRMLIDsCheckBox vertically
  QVBoxLayout * mainLayout = new QVBoxLayout(this);

  // Layout TreeView and (control buttons) horizontally
  QHBoxLayout * horizontalLayout = new QHBoxLayout();
  horizontalLayout->setSpacing(0);
  horizontalLayout->setContentsMargins(0, 0, 0, 0);

  d->TreeView = new QTreeView(this);
  d->TreeView->setModel(d->TreeModel);
  d->TreeView->setRootIsDecorated(false);
  horizontalLayout->addWidget(d->TreeView);

  // Layout control buttons vertically
  QVBoxLayout * controlButtonsLayout = new QVBoxLayout();

  d->ExpandAllButton = new QToolButton(this);
  d->ExpandAllButton->setIcon(QIcon(":/Icons/TreeOpen.png"));
  controlButtonsLayout->addWidget(d->ExpandAllButton);

  d->CollapseAllButton = new QToolButton(this);
  d->CollapseAllButton->setIcon(QIcon(":/Icons/TreeClose.png"));
  controlButtonsLayout->addWidget(d->CollapseAllButton);

  controlButtonsLayout->addItem(
      new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding));

  // Add control buttons to the horizontal layout
  horizontalLayout->addLayout(controlButtonsLayout);

  // Add (TreeView and (control buttons)) to the layout
  mainLayout->addLayout(horizontalLayout);

  d->DisplayMRMLIDsCheckBox = new QCheckBox("Display MRML ID's", this);
  mainLayout->addWidget(d->DisplayMRMLIDsCheckBox);

  // Connect Display MRML Id checkbox
  connect(d->DisplayMRMLIDsCheckBox, SIGNAL(toggled(bool)), SLOT(setMRMLIDsColumnVisible(bool)));

  // Connect control buttons
  connect(d->CollapseAllButton, SIGNAL(clicked()), SLOT(collapseToDepthZero()));
  connect(d->ExpandAllButton, SIGNAL(clicked()), d->TreeView, SLOT(expandAll()));

  // Connect TreeModel
  connect(d->TreeModel, SIGNAL(itemChanged(QStandardItem*)),
          SLOT(onTreeItemChanged(QStandardItem*)));

  // Connect TreeView
//  connect(d->TreeView, SIGNAL(activated(QModelIndex)),
//          SLOT(onTreeItemSelected(QModelIndex)));
  connect(d->TreeView, SIGNAL(clicked(QModelIndex)),
          SLOT(onTreeItemSelected(QModelIndex)));

  this->setStructureNameEditable(false);
  this->setMRMLIDsColumnVisible(false);
  this->setLabelColumnVisible(false);
  this->setClassWeightColumnVisible(false);
  this->setAtlasWeightColumnVisible(false);
  this->setAlphaColumnVisible(false);
}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setMRMLScene(vtkMRMLScene *newScene)
{

}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setup()
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  Q_ASSERT(this->mrmlManager());
  this->populateTreeModel(this->mrmlManager()->GetTreeRootNodeID(),
                          d->TreeModel->invisibleRootItem());

  d->TreeView->expandAll();
}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::collapseToDepthZero()
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  d->TreeView->setUpdatesEnabled(false);
  d->TreeView->collapseAll();
  d->TreeView->expandToDepth(0);
  d->TreeView->setUpdatesEnabled(true);
}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::updateWidgetFromMRML()
{

}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::onTreeItemChanged(QStandardItem * treeItem)
{
  Q_ASSERT(treeItem);
  logger.debug(QString("onTreeItemChanged - %1").arg(treeItem->text()));
  int treeItemType = treeItem->data(ctkPimpl::TreeItemTypeRole).toInt();
  int treeNodeId = treeItem->data(ctkPimpl::TreeNodeIDRole).toInt();

  if (treeItemType == ctkPimpl::StructureNameItemType)
    {
    this->mrmlManager()->SetTreeNodeName(treeNodeId, treeItem->text().toLatin1());
    }
//  else if (treeItem == ctkPimpl::MRMLIDItemType)
//    {

//    }
}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::onTreeItemSelected(const QModelIndex & index)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  QStandardItem * item = d->TreeModel->itemFromIndex(index);
  Q_ASSERT(item);

  int treeNodeId = item->data(ctkPimpl::TreeNodeIDRole).toInt();

  // Get a reference to the associated treeNode
  vtkMRMLEMSTreeNode * currentTreeNode = this->mrmlManager()->GetTreeNode(treeNodeId);
  Q_ASSERT(currentTreeNode);
  if (!currentTreeNode)
    {
    logger.error(QString("onTreeItemSelected - No treeNode associated with id: %1").arg(treeNodeId));
    return;
    }

  emit this->currentTreeNodeChanged(currentTreeNode);

  vtkIdType volumeId = this->mrmlManager()->GetTreeNodeSpatialPriorVolumeID(treeNodeId);
  vtkMRMLVolumeNode * volumeNode = this->mrmlManager()->GetVolumeNode(volumeId);
  emit this->currentSpatialPriorVolumeNodeChanged(volumeNode);
  emit this->currentSpatialPriorVolumeNodeChanged(volumeNode != 0);
}

//-----------------------------------------------------------------------------
namespace
{
void setStructureNameEditableRecursively(QStandardItem * item, bool editable)
{
  Q_ASSERT(item);
  item->setEditable(editable);
  for(int i = 0; i < item->rowCount(); ++i)
    {
    setStructureNameEditableRecursively(item->child(i), editable);
    }
}
}

//-----------------------------------------------------------------------------
CTK_GET_CXX(qSlicerEMSegmentAnatomicalTreeWidget, bool,
            structureNameEditable, StructureNameEditable);

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setStructureNameEditable(bool editable)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  if (d->StructureNameEditable == editable)
    {
    return;
    }

  setStructureNameEditableRecursively(d->TreeModel->invisibleRootItem(), editable);

  d->StructureNameEditable = editable;
}

//-----------------------------------------------------------------------------
bool qSlicerEMSegmentAnatomicalTreeWidget::mrmlIDsColumnVisible() const
{
  CTK_D(const qSlicerEMSegmentAnatomicalTreeWidget);
  return (d->DisplayMRMLIDsCheckBox->checkState() == Qt::Checked);
}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setMRMLIDsColumnVisible(bool visible)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  d->TreeView->header()->setSectionHidden(1, !visible);
  d->DisplayMRMLIDsCheckBox->setChecked(visible);
}

//-----------------------------------------------------------------------------
CTK_GET_CXX(qSlicerEMSegmentAnatomicalTreeWidget, bool, labelColumnVisible, LabelColumnVisible);

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setLabelColumnVisible(bool visible)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  d->TreeView->header()->setSectionHidden(2, !visible);
  d->LabelColumnVisible = visible;
}

//-----------------------------------------------------------------------------
CTK_GET_CXX(qSlicerEMSegmentAnatomicalTreeWidget, bool,
            classWeightColumnVisible, ClassWeightColumnVisible);

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setClassWeightColumnVisible(bool visible)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  d->TreeView->header()->setSectionHidden(3, !visible);
  d->ClassWeightColumnVisible = visible;
}

//-----------------------------------------------------------------------------
CTK_GET_CXX(qSlicerEMSegmentAnatomicalTreeWidget, bool,
            atlasWeightColumnVisible, AtlasWeightColumnVisible);

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setAtlasWeightColumnVisible(bool visible)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  d->TreeView->header()->setSectionHidden(4, !visible);
  d->AtlasWeightColumnVisible = visible;
}

//-----------------------------------------------------------------------------
CTK_GET_CXX(qSlicerEMSegmentAnatomicalTreeWidget, bool, alphaColumnVisible, AlphaColumnVisible);

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::setAlphaColumnVisible(bool visible)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  d->TreeView->header()->setSectionHidden(5, !visible);
  d->AlphaColumnVisible = visible;
}

//-----------------------------------------------------------------------------
void qSlicerEMSegmentAnatomicalTreeWidget::populateTreeModel(
    vtkIdType treeNodeId, QStandardItem * item)
{
  CTK_D(qSlicerEMSegmentAnatomicalTreeWidget);
  Q_ASSERT(this->mrmlManager());
  Q_ASSERT(item);

  // Return if no valid treeNodeId is given
  if (treeNodeId == 0)
    {
    return;
    }

  // Get a reference to the associated treeNode
  vtkMRMLEMSTreeNode * treeNode = this->mrmlManager()->GetTreeNode(treeNodeId);
  Q_ASSERT(treeNode);
  if (!treeNode)
    {
    logger.error(QString("populateTreeModel - No treeNode associated with id: %1").arg(treeNodeId));
    return;
    }
  Q_ASSERT(treeNode->GetParametersNode());

  logger.debug(QString("populateTreeModel - treeNodeId:%1, treeNodeName: %2").
               arg(treeNodeId).arg(treeNode->GetName()));

  QList<QStandardItem*> itemList;

  // Structure item
  QStandardItem * structureItem = new QStandardItem(QString("%1").arg(treeNode->GetName()));
  structureItem->setData(QVariant(treeNodeId), ctkPimpl::TreeNodeIDRole);
  structureItem->setData(QVariant(ctkPimpl::StructureNameItemType), ctkPimpl::TreeItemTypeRole);
  structureItem->setEditable(d->StructureNameEditable);
  itemList << structureItem;

  // MRML ID item
  QStandardItem * mrmlIDItem = new QStandardItem(QString("%1").arg(treeNode->GetID()));
  mrmlIDItem->setData(QVariant(treeNodeId), ctkPimpl::TreeNodeIDRole);
  mrmlIDItem->setData(QVariant(ctkPimpl::MRMLIDItemType), ctkPimpl::TreeItemTypeRole);
  mrmlIDItem->setEditable(false);
  itemList << mrmlIDItem;

  // Label item - Available only for tree leaf
  QStandardItem * labelItem = new QStandardItem();
  labelItem->setData(QVariant(treeNodeId), ctkPimpl::TreeNodeIDRole);
  labelItem->setEditable(false);
  if (treeNode->GetNumberOfChildNodes() == 0) // Is treeNode a leaf ?
    {
    Q_ASSERT(treeNode->GetParametersNode()->GetLeafParametersNode());
    // TODO label should be editable
    labelItem->setText(QString("%1").arg(
        treeNode->GetParametersNode()->GetLeafParametersNode()->GetIntensityLabel()));
    labelItem->setData(QVariant(ctkPimpl::LabelItemType), ctkPimpl::TreeItemTypeRole);
    }
  itemList << labelItem;

  // ClassWeight item
  QStandardItem * classWeightItem = new QStandardItem(QString("%1").arg(
      treeNode->GetParametersNode()->GetClassProbability()));
  classWeightItem->setEditable(false); // TODO class weight should be editable
  classWeightItem->setData(QVariant(treeNodeId), ctkPimpl::TreeNodeIDRole);
  classWeightItem->setData(QVariant(ctkPimpl::ClassWeightItemType), ctkPimpl::TreeItemTypeRole);
  itemList << classWeightItem;

  // AtlasWeight item
  QStandardItem * atlasWeightItem = new QStandardItem(QString("%1").arg(
      treeNode->GetParametersNode()->GetSpatialPriorWeight()));
  atlasWeightItem->setEditable(false); // TODO atlas weight should be editable
  atlasWeightItem->setData(QVariant(treeNodeId), ctkPimpl::TreeNodeIDRole);
  atlasWeightItem->setData(QVariant(ctkPimpl::AtlasWeightItemType), ctkPimpl::TreeItemTypeRole);
  itemList << atlasWeightItem;

  // Alpha item - Available only for none tree leaf
  QStandardItem * alphaItem = new QStandardItem();
  alphaItem->setData(QVariant(treeNodeId), ctkPimpl::TreeNodeIDRole);
  alphaItem->setEditable(false);
  if (treeNode->GetNumberOfChildNodes() != 0) // Is treeNode NOT a leaf ?
    {
    Q_ASSERT(treeNode->GetParametersNode()->GetParentParametersNode());
    alphaItem->setText(QString("%1").arg(
        treeNode->GetParametersNode()->GetParentParametersNode()->GetAlpha()));
    // TODO alpha should be editable
    alphaItem->setData(QVariant(ctkPimpl::AlphaItemType), ctkPimpl::TreeItemTypeRole);
    }
  itemList << alphaItem;


  item->appendRow(itemList);

  // Loop through current node children and recursively call ourself
  int numberOfChildren = this->mrmlManager()->GetTreeNodeNumberOfChildren(treeNodeId);
  for (int i = 0; i < numberOfChildren; i++)
    {
    this->populateTreeModel(
      this->mrmlManager()->GetTreeNodeChildNodeID(treeNodeId, i), structureItem);
    }
}
