#include "tabledelegate.h"

#include <QLineEdit>
#include <QDoubleValidator>

TableDelegate::TableDelegate(QObject *parent) :
    QItemDelegate(parent)
{
}

QWidget *TableDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    QLineEdit *editor = new QLineEdit(parent);
    editor->setValidator(new QDoubleValidator);
    return editor;
}


//void TableDelegate::setEditorData(QWidget *editor,
//                                 const QModelIndex &index) const
//{
//    QString value =index.model()->data(index, Qt::EditRole).toString();
//        QLineEdit *line = static_cast<QLineEdit*>(editor);
//        line->setText(value);
//}


//void TableDelegate::setModelData(QWidget *editor,
//                                QAbstractItemModel *model,
//                                const QModelIndex &index) const
//{
//    QLineEdit *line = static_cast<QLineEdit*>(editor);
//    QString value = line->text();
//    model->setData(index, value);
//}


//void TableDelegate::updateEditorGeometry(QWidget *editor,
//                                        const QStyleOptionViewItem &option,
//                                        const QModelIndex &index) const
//{
//    editor->setGeometry(option.rect);
//}
