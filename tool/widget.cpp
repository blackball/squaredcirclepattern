#include "widget.h"
#include "ui_widget.h"
#include <QPainter>
#include <QBrush>
#include <QFileDialog>

Widget::Widget(QWidget *parent) :
        QWidget(parent),
        ui(new Ui::Widget)
{
        ui->setupUi(this);
}

Widget::~Widget()
{
        delete ui;
}

void Widget::on_pushButton_clicked()
{
        const QString type = ui->comboBox->currentText();
        if (type == QString("Ours")) {
                drawOursPattern();
        }
        else if (type == QString("ChessBoard")) {
                drawChessBoardPattern();
        }
        else if(type == QString("Circle")) {
                drawCirclePattern();
        }
        else if (type == QString("Ring")) {
                drawRingPattern();
        }
}

void Widget::drawOursPattern()
{
        int rows = ui->rowsLineEdit->text().toInt();
        int cols = ui->colsLineEdit->text().toInt();
        int size = ui->sizeLineEdit->text().toInt();

        const int ox = 1;
        const int oy = 1;

        int radius = size/4;

        int width = cols * size + (cols+1) * radius + size;
        int height = rows * size + (rows+1) * radius + size;

        QImage img(width, height, QImage::Format_RGB888);
        img.fill(Qt::white);

        drawRect(img, QRect(radius, radius, width - 2*radius, height- 2*radius), Qt::black);
        drawRect(img, QRect(2*radius, 2*radius, width - 4*radius, height - 4*radius), Qt::white);


        for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                        if (i != ox || j != oy) {
                                QRect r(2*radius + radius * (j+1) + j * size, 2*radius + radius * (i+1) + i*size, size, size);
                                QPoint center = r.center();
                                drawRect(img, r, Qt::black);
                                drawCircle(img, center, radius, Qt::white);
                                //drawCircle(img, center, 3, Qt::black);
                        }
                        else {
                                QRect r(2*radius + radius * (j+1) + j * size, 2*radius + radius * (i+1) + i*size, size, size);
                                QPoint center = r.center();
                                drawRect(img, r, Qt::black);
                                drawCircle(img, center, radius*1.5, Qt::white);
                                //drawCircle(img, center, 3, Qt::black);
                        }
                }
        }

        m_image = img;
        ui->imageLabel->setPixmap(QPixmap::fromImage(img.scaled(ui->imageLabel->size())));
        this->update();
}

void Widget::drawChessBoardPattern()
{
        int rows = ui->rowsLineEdit->text().toInt();
        int cols = ui->colsLineEdit->text().toInt();
        int size = ui->sizeLineEdit->text().toInt();

        int width = (cols + 2) * size;
        int height = (rows + 2) * size;

        QImage img(width, height, QImage::Format_RGB888);
        img.fill(Qt::white);

        bool alter = false;
        for (int i = 1; i < rows+1; ++i) {
                for (int j = 1; j < cols+1; ++j) {
                        if (alter) {
                                QRect r(j*size, i*size, size, size);
                                drawRect(img, r, Qt::black);
                        }
                        alter = !alter;
                }
                alter = !alter;
        }

        m_image = img;

        ui->imageLabel->setPixmap(QPixmap::fromImage(img.scaled(ui->imageLabel->size())));
        this->update();
}

void Widget::drawCirclePattern()
{
        // TODO
}

void Widget::drawRingPattern()
{
        // TODO
}

void Widget::drawRect(QImage &img, QRect rect, QColor fillColor)
{
        QPainter painter(&img);
        painter.setRenderHint(QPainter::Antialiasing);
        painter.setPen(Qt::white);
        painter.setBrush(QBrush(fillColor));
        painter.drawRect(rect);
}

void Widget::drawCircle(QImage &img, QPoint origin, int radius, QColor fillColor)
{
        QPainter painter(&img);
        painter.setRenderHint(QPainter::Antialiasing);
        painter.setPen(Qt::white);
        painter.setBrush(QBrush(fillColor));
        painter.drawEllipse(origin, radius, radius);
}

void Widget::on_pushButton_2_clicked()
{
        QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                        "untitled.png",
                                                        tr("Images (*.png *.bmp *.jpg)"));
        if (fileName.isNull()) {
                return ;
        }

        if (m_image.isNull()) {
                return ;
        }

        m_image.save(fileName);
}
