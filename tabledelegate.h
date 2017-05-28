#ifndef TABLEDELEGATE_H
#define TABLEDELEGATE_H

#include <QItemDelegate>

class TableDelegate : public QItemDelegate
{
    Q_OBJECT
    public:
        explicit TableDelegate(QObject *parent = 0);

    protected:
        QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
//        void setEditorData(QWidget * editor, const QModelIndex & index) const;
//        void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const;
//        void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const;

};

#endif // TABLEDELEGATE_H
