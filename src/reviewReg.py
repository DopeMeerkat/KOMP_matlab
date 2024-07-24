import sys
import os
import itertools
from PyQt5 import QtCore, QtGui, QtWidgets
from PIL import Image

# TODO: 
# - loaded image != viewing image
# - path before directory contains '_type'
# - iterate backwards
# - left and right
# - resolution decreased?

Image.MAX_IMAGE_PIXELS = 200000000
IMAGE_HEIGHT = 1100 # minimum height of image in pixels, change if needed based on device
IMAGE_WIDTH = 600  # minimum width of image in pixels

class GraphicView(QtWidgets.QGraphicsView):
    rectChanged = QtCore.pyqtSignal(QtCore.QRect)

    def __init__(self, *args, **kwargs):
        QtWidgets.QGraphicsView.__init__(self, *args, **kwargs)
        self.cursor = None
        # self.cursor = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, self)
        self.setMouseTracking(True)
        # self.origin = QtCore.QPoint()
        # self.changeRubberBand = False
        self.scene = QtWidgets.QGraphicsScene()
        self.setScene(self.scene)
        self.scene.setSceneRect(QtCore.QRectF(0,0,IMAGE_WIDTH, IMAGE_HEIGHT))
        # self.selectedRegion = {'x':0,'y':0,'w':-1,'h':-1}
        self.graphicsPixmapItem = None

        self.drawing = False
        self.lastPoint = QtCore.QPoint()
        self.penColor = QtCore.Qt.white
        self.penSize = 18
        self.pen = None
        self.brush = None
        self.updatePen()

        
    def updatePen(self):
        self.pen = QtGui.QPen(self.penColor, 1, QtCore.Qt.SolidLine)
        self.brush = QtGui.QBrush(self.penColor)
        # if self.penColor == QtCore.Qt.black:
        #     self.cursor = QtGui.QCursor(QtGui.QPixmap('blackCursor.png').scaled(int(self.penSize/2), int(self.penSize/2)))
        # else:
        #     self.cursor = QtGui.QCursor(QtGui.QPixmap('whiteCursor.png').scaled(int(self.penSize/2), int(self.penSize/2)))
        self.cursor = QtGui.QCursor(QtGui.QPixmap('whiteCursor.png').scaled(int(self.penSize/2), int(self.penSize/2)))
        self.setCursor(self.cursor)
        # self.cursor.setGeometry(QtCore.QRect(0,0, self.penSize, self.penSize))

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.drawing = True
            self.lastPoint = event.pos()

    def mouseMoveEvent(self, event):
        if event.buttons() and QtCore.Qt.LeftButton and self.drawing:
            # self.scene.addLine(QtCore.QLineF(self.lastPoint, event.pos()), self.pen)
            self.scene.addEllipse(event.pos().x() - int(self.penSize / 2), 
                                event.pos().y() - int(self.penSize / 2), 
                                self.penSize, self.penSize, self.pen, brush = self.brush)
            self.lastPoint = event.pos()
            self.update()
        
        # self.cursor.setGeometry(QtCore.QRect(event.pos().x() - int(self.penSize / 2), 
        #                                         event.pos().y() - int(self.penSize / 2), 
        #                                         self.penSize, self.penSize))

    def mouseReleaseEvent(self, event):
        if event.button == QtCore.Qt.LeftButton:
            self.drawing = False

    # def mousePressEvent(self, event):
    #     self.origin = event.pos()
    #     self.rubberBand.setGeometry(QtCore.QRect(self.origin, QtCore.QSize()))
    #     self.rectChanged.emit(self.rubberBand.geometry())
    #     self.rubberBand.show()
    #     self.changeRubberBand = True
    #     QtWidgets.QGraphicsView.mousePressEvent(self, event)

    # def mouseMoveEvent(self, event):
    #     if self.changeRubberBand:
    #         self.rubberBand.setGeometry(QtCore.QRect(self.origin, event.pos()).normalized())
    #         self.rectChanged.emit(self.rubberBand.geometry())
    #     QtWidgets.QGraphicsView.mouseMoveEvent(self, event)

    # def mouseReleaseEvent(self, event):
    #     self.changeRubberBand = False
    #     QtWidgets.QGraphicsView.mouseReleaseEvent(self, event)
    #     self.selectedRegion['x'] = self.rubberBand.geometry().x()
    #     self.selectedRegion['y'] = self.rubberBand.geometry().y()
    #     self.selectedRegion['w'] = self.rubberBand.geometry().width()
    #     self.selectedRegion['h'] = self.rubberBand.geometry().height()
    #     self.scene.addRect(self.selectedRegion['x'],self.selectedRegion['y'],self.selectedRegion['w'],self.selectedRegion['h'], pen = QtGui.QPen(QtCore.Qt.red, 4))
    #     self.rubberBand.hide()

class ImageLoader(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        layout = QtWidgets.QGridLayout(self)

        self.setWindowTitle('Cropping Tool')
        self.label = GraphicView()
        layout.addWidget(self.label, 1, 2, 1, 2)
        self.label.setMinimumSize(IMAGE_WIDTH, IMAGE_HEIGHT)
        self.label.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.label.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.label.setAlignment(QtCore.Qt.AlignCenter)

        self.label2 = QtWidgets.QGraphicsView()
        self.label2.scene = QtWidgets.QGraphicsScene()
        self.label2.setScene(self.label2.scene)
        self.label2.scene.setSceneRect(QtCore.QRectF(0,0,IMAGE_WIDTH, IMAGE_HEIGHT))
        self.label2.graphicsPixmapItem = None
        layout.addWidget(self.label2, 1, 4, 1, 2)
        self.label2.setMinimumSize(IMAGE_WIDTH, IMAGE_HEIGHT)
        self.label2.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.label2.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.label2.setAlignment(QtCore.Qt.AlignCenter)

        self.loadImageButton = QtWidgets.QPushButton('Load image')
        layout.addWidget(self.loadImageButton, 0, 0, 1, 2)

        self.clearButton = QtWidgets.QPushButton('Clear')
        layout.addWidget(self.clearButton, 0, 2, 1, 1)

        self.scrapButton = QtWidgets.QPushButton('Scrap sample')
        layout.addWidget(self.scrapButton, 0, 3, 1, 1)

        self.penWhiteButton = QtWidgets.QPushButton('White')
        layout.addWidget(self.penWhiteButton, 2, 2, 1, 1)

        self.penBlackButton = QtWidgets.QPushButton('Black')
        layout.addWidget(self.penBlackButton, 2, 3, 1, 1)
        
        self.penUpButton = QtWidgets.QPushButton('Bigger')
        layout.addWidget(self.penUpButton, 2, 4, 1, 1)

        self.penDownButton = QtWidgets.QPushButton('Smaller')
        layout.addWidget(self.penDownButton, 2, 5, 1, 1)

        self.saveChangesButton = QtWidgets.QPushButton('Save changes')
        layout.addWidget(self.saveChangesButton, 0, 6, 1, 2)



        # self.actionLabel = QtWidgets.QLabel()
        # layout.addWidget(self.actionLabel, 2, 2, 1, 2)
        # self.actionLabel.setText("Load a file to get started")
        # self.actionLabel.setAlignment(QtCore.Qt.AlignCenter)

        self.nextImageButton = QtWidgets.QPushButton('>')
        layout.addWidget(self.nextImageButton, 1, 6)
        self.nextImageButton.setMinimumSize(20,IMAGE_HEIGHT)

        self.prevImageButton = QtWidgets.QPushButton('<')
        layout.addWidget(self.prevImageButton, 1, 1)
        self.prevImageButton.setMinimumSize(20,IMAGE_HEIGHT)


        self.loadImageButton.clicked.connect(self.loadImage)
        self.clearButton.clicked.connect(self.clearScene)
        self.scrapButton.clicked.connect(self.scrap)
        self.penWhiteButton.clicked.connect(self.penWhite)
        self.penBlackButton.clicked.connect(self.penBlack)
        self.penUpButton.clicked.connect(self.penUp)
        self.penDownButton.clicked.connect(self.penDown)
        self.saveChangesButton.clicked.connect(self.cropImage)
        self.nextImageButton.clicked.connect(self.nextImage)
        self.prevImageButton.clicked.connect(self.prevImage)

        self.filename = ''
        self.filename2 = ''
        self.dirIterator = None
        self.regIterator = None
        self.fileList = []
        self.regList = []
        self.pixmap = None
        self.pixmap2 = None
        self.dirname = ''

    def loadImage(self):
        self.clearScene()
        self.dirname = os.path.dirname(os.getcwd())
        # self.dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Folder')
        self.fileList = []
        self.regListList = []
        layerPath = os.path.join(self.dirname, os.path.join('01_Submitted', 'Layers'))
        thresholdPath = os.path.join(layerPath,'Threshold')
        regPath = os.path.join(layerPath,'Images')

        # self.setWindowTitle(self.filename)
        # self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
        #     QtCore.Qt.KeepAspectRatio)
        # if self.pixmap.isNull():
        #     return
        # # self.label.graphicsPixmapItem = QtWidgets.QGraphicsPixmapItem(QtGui.QPixmap(self.pixmap))
        # # self.label.scene.addItem(self.label.graphicsPixmapItem)
        # self.clearScene()
        # dirpath = os.path.dirname(self.filename)
        
        for f in os.listdir(thresholdPath):
            fpath = os.path.join(thresholdPath, f)
            if f.find('_2_sh') != -1:
                self.fileList.append(fpath)

        for f in os.listdir(regPath):
            fpath = os.path.join(regPath, f)
            if f.find('2_NoDAPI') != -1:
                self.regList.append(fpath)


        self.fileList.sort()
        self.regList.sort()
        self.dirIterator = iter(self.fileList)
        self.regIterator = iter(self.regList)
        self.filename = next(self.dirIterator)
        self.filename2 = next(self.regIterator)

        self.setWindowTitle(os.path.basename(self.filename) + ' --- ' + os.path.basename(self.filename2))

        # print('fileList', len(self.fileList))
        # print('regList', len(self.regList))
        # for f in self.fileList:
        #     print(f)
        # for f in self.regList:
        #     print(f)
        self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
            QtCore.Qt.KeepAspectRatio)
        self.pixmap2 = QtGui.QPixmap(self.filename2).scaled(self.label2.size(), 
            QtCore.Qt.KeepAspectRatio)
        if self.pixmap.isNull() or self.pixmap2.isNull():
            return
        self.clearScene()

        # try:
        #     while True:
        #         # cycle through the iterator until the current file is found
        #         if os.path.basename(next(self.dirIterator)) == os.path.basename(self.filename): # windows uses '\'
        #             break
        # except:
        #     # self.filename = next(self.dirIterator)
        #     self.actionLabel.setText('Warning: Channel being viewed is not the designated one for the type')

        

    def nextImage(self):
        if self.fileList:
            try:
                self.filename = next(self.dirIterator)
                self.filename2 = next(self.regIterator)
                self.setWindowTitle(os.path.basename(self.filename) + ' --- ' + os.path.basename(self.filename2))
                self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
                    QtCore.Qt.KeepAspectRatio)
                self.pixmap2 = QtGui.QPixmap(self.filename2).scaled(self.label2.size(), 
                    QtCore.Qt.KeepAspectRatio)
                if self.pixmap.isNull():
                    # the file is not a valid image, remove it from the list
                    # and try to load the next one
                    self.fileList.remove(self.filename)
                    self.nextImage()
                else:
                    self.clearScene()
            except:
                # the iterator has finished, restart it
                self.dirIterator = iter(self.fileList)
                self.regIterator = iter(self.regList)
                self.nextImage()
        else:
            # no file list found, load an image
            self.loadImage()

    def prevImage(self):
        if self.fileList:
            for _ in itertools.repeat(None, len(self.fileList) - 1):
                try:
                    self.filename = next(self.dirIterator)
                    self.filename2 = next(self.regIterator)
                    
                except:
                    # the iterator has finished, restart it
                    self.dirIterator = iter(self.fileList)
                    self.regIterator = iter(self.regList)
                    self.nextImage()

            self.setWindowTitle(os.path.basename(self.filename) + ' --- ' + os.path.basename(self.filename2))
            self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
                QtCore.Qt.KeepAspectRatio)
            self.pixmap2 = QtGui.QPixmap(self.filename2).scaled(self.label2.size(), 
                QtCore.Qt.KeepAspectRatio)
            # if self.pixmap.isNull():
            #     # the file is not a valid image, remove it from the list
            #     # and try to load the next one
            #     self.fileList.remove(self.filename)
            #     self.nextImage()
            # else:
            self.clearScene()
        else:
            # no file list found, load an image
            self.loadImage()


    def clearScene(self):
        # print('clearing')
        self.label.scene.clear()
        self.label.setAlignment(QtCore.Qt.AlignTop)
        self.label.setAlignment(QtCore.Qt.AlignLeft)
        self.label.graphicsPixmapItem = QtWidgets.QGraphicsPixmapItem(QtGui.QPixmap(self.pixmap))
        self.label.scene.addItem(self.label.graphicsPixmapItem)

        self.label2.scene.clear()
        self.label2.setAlignment(QtCore.Qt.AlignTop)
        self.label2.setAlignment(QtCore.Qt.AlignLeft)
        self.label2.graphicsPixmapItem = QtWidgets.QGraphicsPixmapItem(QtGui.QPixmap(self.pixmap2))
        self.label2.scene.addItem(self.label2.graphicsPixmapItem)

    def penBlack(self):
        self.label.penColor = QtCore.Qt.black
        self.label.updatePen()
    
    def penWhite(self):
        self.label.penColor = QtCore.Qt.white
        self.label.updatePen()

    def penUp(self):
        self.label.penSize += 2
        self.label.updatePen()

    def penDown(self):
        self.label.penSize -= 2
        self.label.updatePen()

    def scrap(self):
        if os.path.exists(os.path.join(self.dirname,'scraplist.txt')):
            file = open(os.path.join(self.dirname,'scraplist.txt'), 'a')
        else:
            file = open(os.path.join(self.dirname,'scraplist.txt'), 'w+')
        basename = os.path.basename(self.filename)
        i = basename.find('_hF_')
        g = basename[i + 4]
        l = basename[i + 6]
        s = basename[i + 9]
        file.write(f'{g} {l} {s}\n')
        file.close()

    def cropImage(self):
        areaF = self.label.scene.sceneRect()
        image = QtGui.QImage(self.filename)
        # image = QtGui.QImage(areaF.toRect().size(),QtGui.QImage.Format_ARGB32_Premultiplied)
        # dpm = int(300 / 0.0254)
        # image.setDotsPerMeterY(dpm)
        # image.setDotsPerMeterX(dpm)
        painter = QtGui.QPainter(image)
        QtCore.QRectF(image.rect())

        self.label.scene.render(painter, QtCore.QRectF(image.rect()), areaF, QtCore.Qt.KeepAspectRatioByExpanding)
        # print(QtCore.QRectF(image.rect()))
        # print(areaF)
        painter.end()
        # image.save(os.path.join(self.dirname, 'test.jpg'), format='JPEG', quality=100)
        if not os.path.exists(self.filename[:-4] + '_original.jpg'):
            os.rename(self.filename, self.filename[:-4] + '_original.jpg')
        else:
            os.remove(self.filename)
        image.save(self.filename, format='JPEG', quality=100)
        self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
                QtCore.Qt.KeepAspectRatio)
        self.clearScene()

        # base = self.filename[:self.filename.find('_s') + 3]
        # for f in self.fullFileList:
        #     if f.find(base) != -1: #contains base
        #         # print('base =', base, 'cropping', f)
        #         im = Image.open(f)
        #         # print('im format:', im.format, 'dpi', im.info['dpi'])
        #         width, height = im.size
        #         ratio = width / self.pixmap.width()
        #         x1 = int(self.label.selectedRegion['x'] * ratio)
        #         y1 = int(self.label.selectedRegion['y'] * ratio)
        #         x2 = min(width, int((self.label.selectedRegion['w'] + self.label.selectedRegion['x']) * ratio))
        #         y2 = min(height, int((self.label.selectedRegion['h'] + self.label.selectedRegion['y']) * ratio))
        #         im1 = im.crop((x1,y1,x2,y2))
        #         # im1 = im.crop((int(self.label.selectedRegion['x'] * ratio), 
        #         #                int(self.label.selectedRegion['y'] * ratio), 
        #         #                int((self.label.selectedRegion['w'] + self.label.selectedRegion['x']) * ratio), 
        #         #                int((self.label.selectedRegion['h'] + self.label.selectedRegion['y']) * ratio)))
        #         # print('im1 format:', im.format, im.info['dpi'])
        #         im1.save(f, format = 'JPEG', dpi = im1.info['dpi'])#, quality = 'keep')
        # self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
        #             QtCore.Qt.KeepAspectRatio)
        # self.clearScene()

        

        



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    imageLoader = ImageLoader()
    imageLoader.show()
    sys.exit(app.exec_())