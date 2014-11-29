#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

namespace Ui {
        class Widget;
}

class Widget : public QWidget
{
        Q_OBJECT
public:
        explicit Widget(QWidget *parent = 0);
        ~Widget();
                 
private slots:
        void on_pushButton_clicked();
        void on_pushButton_2_clicked();

private:
        void drawOursPattern();
        void drawChessBoardPattern();
        void drawCirclePattern();
        void drawRingPattern();
        
        void drawRect(QImage &img, QRect rect, QColor fillColor);
        void drawCircle(QImage &img, QPoint origin, int radius, QColor fillColor);

private:
        Ui::Widget *ui;
        QImage m_image;
};

#endif // WIDGET_H
