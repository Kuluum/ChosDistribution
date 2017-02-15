#ifndef FILEREADER_H
#define FILEREADER_H

#include <QObject>

class FileReader : public QObject
{
    Q_OBJECT

public:
    explicit FileReader(QObject *parent = 0);
    QStringList readFileWithDialog(QWidget *parentWidget);
    QStringList readFileNamed(QString fileName);

private:

signals:

public slots:
};

#endif // FILEREADER_H
