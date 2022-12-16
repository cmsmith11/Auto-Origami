#!/usr/bin/env python3

# from PIL import Image

# print('Starting stickman.py')

# im = Image.open('BW_man.png')
# w, h = im.size
# print(im.size)
# print(w, h)
# px = im.load()
# print(px[w-1, h-1])


# def lum(pixel):
#   return (pixel[0] + pixel[1] + pixel[2]) / 3

# seen = set()
# values = []

# for r in range(im.size[1]):
#   values.append([])
#   for c in range(im.size[0]):
#       values[-1].append(0)


# for r in range(im.size[1]):
#   for c in range(im.size[0]):
#       # loop over every pixel
#       rec()


# px[w//2, h//2] = (255, 255, 255, 255)
# print(px[w//2, h//2])



# im.show()


# import tkinter as tk

# print('starting stickman')

# clickCount = 0 # used for checking start or end of line
# points = []

# def click(event):
#   global clickCount
#   global points

#   clickCount += 1
#   point = Point(event.x, event.y)

#   if clickCount % 2 != 0:
#       # check if close to an existing point. If so, assign that to point
#       pass
#   else:
#       # Make connection from new point to previous point
#       pass

#   # draw stuff

#   print(event)

# class Point:
#   def __init__(self, x, y):
#       self.x = x
#       self.y = y
#       self.neighbors = [] # list of Points which are connected to this point

# root = tk.Tk()
# frame = tk.Frame(root, width=600, height=600)
# frame.config(bg='white')
# frame.bind("<Button-1>", click)
# frame.pack()
# root.mainloop()




# import sys
# from PyQt5.QtCore import *
# from PyQt5.QtGui import *
# from PyQt5.QtWidgets import *

# clickCount = 0 # used for checking start or end of line
# points = []

# class ClickWidget(QWidget):
#     def __init__(self, mw):
#         super().__init__()
#         self.mw = mw

#     def mousePressEvent(self, event):
#         print(event.pos().x(), event.pos().y())
#         self.mw.draw_something()


# def main():
#    app = QApplication(sys.argv)
#    w = ClickWidget()
#    w.setGeometry(100,100,600,600)
#    w.setStyleSheet("background-color: white;")
#    w.setWindowTitle("Make Stick Figure")
#    print('start')
#    w.show()
#    print('end?')
#    sys.exit(app.exec_()) # mainloop

# if __name__ == '__main__':
#    main()


import sys
import math
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import Qt

# class ClickWidget(QtWidgets.QWidget):
#     def __init__(self, mw):
#         super().__init__()
#         self.mw = mw

#     def mousePressEvent(self, event):
#         print(event.pos().x(), event.pos().y())
#         self.mw.draw_something()

def writeToFile(window, filename):
    f = open(filename, "w")
    f.write(window.__str__())
    f.close()

class Point:
    def __init__(self, x, y):
        self.i = -1
        self.x = x
        self.y = y
        self.neighbors = [] # list of Points which are connected to this point
    def dist(self, p):
        return math.sqrt((p.x - self.x)**2 + (p.y - self.y)**2)
    def __str__(self):
        return('(' + str(self.x) + ', ' + str(self.y) + ')')


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        canvas = QtGui.QPixmap(600, 600)

        self.points = []
        self.lastPoint = None # technically 2nd to last, when used
        self.clickCount = 0
        self.label = QtWidgets.QLabel()

        self.label.setPixmap(canvas)
        self.setCentralWidget(self.label)

        self.b1 = QtWidgets.QPushButton("Done", self)
        print(canvas.rect().width(), canvas.rect().height())
        self.w = canvas.rect().width()
        self.h = canvas.rect().height()
        self.b1.resize(100, 32)
        self.b1.move(self.w - 125, self.h - 50) 
        self.b1.clicked.connect(self.onClick)

    # format: x y n1,n2,n3
    def __str__(self):
        # assign all points an index
        for i, p in enumerate(self.points):
            p.i = i

        s = ''
        for p in self.points:
            s += str(p.x / self.w) + ' ' + str((self.h - p.y) / self.h) + ' '
            for n in p.neighbors:
                s += str(n.i) + ','
            s = s[:-1] # remove last ','
            s += '\n'
        return(s)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Enter - 1:
            self.deleteLater()

    def onClick(self):
        self.deleteLater()

    def mousePressEvent(self, event):
        self.clickCount += 1
        point = Point(event.pos().x(), event.pos().y())
        if self.clickCount % 2 != 0:
            self.lastPoint = point
            # check if close to another point
            nearest = None
            for p in self.points:
                if point.dist(p) < 30:
                    print('same point!!')
                    if nearest == None or point.dist(p) < point.dist(nearest):
                        nearest = p
            if nearest != None:
                point = nearest
                self.lastPoint = point
            else:
                self.points.append(point)
        else:
            self.points.append(point)
        self.draw()

    def draw(self):
        painter = QtGui.QPainter(self.label.pixmap())
        pen = QtGui.QPen()
        pen.setWidth(4)
        painter.setPen(pen)
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        ps = self.points
        p1 = ps[-1]
        painter.drawPoint(ps[-1].x, ps[-1].y)
        if self.clickCount % 2 == 0:
            p2 = self.lastPoint
            p1.neighbors.append(p2)
            p2.neighbors.append(p1)
            painter.drawLine(p1.x, p1.y, p2.x, p2.y)
        self.update()

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()
    writeToFile(window, "temp.stk")

if __name__ == '__main__':
    main()
    print('end?')
