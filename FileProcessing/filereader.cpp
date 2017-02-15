#include "FileProcessing/filereader.h"
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>
FileReader::FileReader(QObject *parent) : QObject(parent)
{

}

QStringList FileReader::readFileWithDialog(QWidget *ui) {

   QString fileName = QFileDialog::getOpenFileName(ui,
                                 tr("Open Document"),
                                 QDir::currentPath(),
                                 tr("Text files (*.txt);;All files (*.*)"),
                                 0,
                                 QFileDialog::DontUseNativeDialog );

       return readFileNamed(fileName);
}

QStringList FileReader::readFileNamed(QString fileName) {
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(0, "error", file.errorString());
    }

    QTextStream in(&file);
    QStringList stringList;

    while(!in.atEnd()) {
        QString line = in.readLine();
        stringList.append(line);
    }

    file.close();

    return stringList;
}
