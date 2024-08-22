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
IMAGE_HEIGHT = 1200 # minimum height of image in pixels, change if needed based on device
IMAGE_WIDTH = 900  # minimum width of image in pixels

class BB():
    def __init__(self, x, y, w, h):
        self.left = x
        self.right = x + w
        self.top = y
        self.bot = y + h

class GraphicView(QtWidgets.QGraphicsView):
    rectChanged = QtCore.pyqtSignal(QtCore.QRect)

    def __init__(self, *args, **kwargs):
        QtWidgets.QGraphicsView.__init__(self, *args, **kwargs)
        self.rubberBand = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, self)
        self.setMouseTracking(True)
        self.origin = QtCore.QPoint()
        self.changeRubberBand = False
        self.scene = QtWidgets.QGraphicsScene()
        self.scene.setSceneRect(QtCore.QRectF(0,0,IMAGE_WIDTH, IMAGE_HEIGHT))
        self.setScene(self.scene)
        self.selectedRegion = {'x':0,'y':0,'w':-1,'h':-1}
        self.BBList = []
        self.graphicsPixmapItem = None
        self.minPixels = 0


    def mousePressEvent(self, event):
        self.origin = event.pos()
        self.rubberBand.setGeometry(QtCore.QRect(self.origin, QtCore.QSize()))
        self.rectChanged.emit(self.rubberBand.geometry())
        self.rubberBand.show()
        self.changeRubberBand = True
        QtWidgets.QGraphicsView.mousePressEvent(self, event)

    def mouseMoveEvent(self, event):
        if self.changeRubberBand:
            self.rubberBand.setGeometry(QtCore.QRect(self.origin, event.pos()).normalized())
            self.rectChanged.emit(self.rubberBand.geometry())
        QtWidgets.QGraphicsView.mouseMoveEvent(self, event)

    def mouseReleaseEvent(self, event):
        self.changeRubberBand = False
        QtWidgets.QGraphicsView.mouseReleaseEvent(self, event)
        self.selectedRegion['x'] = self.rubberBand.geometry().x()
        self.selectedRegion['y'] = self.rubberBand.geometry().y()
        self.selectedRegion['w'] = self.rubberBand.geometry().width()
        self.selectedRegion['h'] = self.rubberBand.geometry().height()
        # print(len(self.BBList))
        if len(self.BBList) < 2:
            self.scene.addRect(self.selectedRegion['x'],self.selectedRegion['y'],self.selectedRegion['w'],self.selectedRegion['h'], pen = QtGui.QPen(QtCore.Qt.red, 4))
            self.BBList.append(BB(self.selectedRegion['x'],self.selectedRegion['y'],self.selectedRegion['w'],self.selectedRegion['h']))
        self.rubberBand.hide()


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

        self.loadImageButton = QtWidgets.QPushButton('Load image')
        layout.addWidget(self.loadImageButton, 0, 0, 1, 2)

        self.clearButton = QtWidgets.QPushButton('Clear')
        layout.addWidget(self.clearButton, 0, 2, 1, 1)

        self.saveChangesButton = QtWidgets.QPushButton('Save changes')
        layout.addWidget(self.saveChangesButton, 0, 4, 1, 2)

        self.nextImageButton = QtWidgets.QPushButton('>')
        layout.addWidget(self.nextImageButton, 1, 4)
        self.nextImageButton.setMinimumSize(20,IMAGE_HEIGHT)

        self.prevImageButton = QtWidgets.QPushButton('<')
        layout.addWidget(self.prevImageButton, 1, 1)
        self.prevImageButton.setMinimumSize(20,IMAGE_HEIGHT)


        self.loadImageButton.clicked.connect(self.loadImage)
        self.clearButton.clicked.connect(self.clearScene)
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
        self.loadImage()

    def loadImage(self):
        self.clearScene()
        self.fileList = []
        self.regListList = []

        # for testing
        self.dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Folder')
        errorPath = self.dirname

        # self.dirname = os.getcwd()
        # layerPath = os.path.join(self.dirname, os.path.join('01_Submitted', 'Layers'))
        # errorPath = os.path.join(layerPath,'Images')
        
        for f in os.listdir(errorPath):
            fpath = os.path.join(errorPath, f)
            if f.find('_error.') != -1:
                self.fileList.append(fpath)

        self.fileList.sort()
        self.regList.sort()
        self.dirIterator = iter(self.fileList)
        self.filename = next(self.dirIterator)

        self.setWindowTitle(os.path.basename(self.filename))

        self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
            QtCore.Qt.KeepAspectRatio)
        if self.pixmap.isNull():
            return
        self.clearScene()

    def nextImage(self):
        if self.fileList:
            try:
                self.filename = next(self.dirIterator)
                self.setWindowTitle(os.path.basename(self.filename))
                self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
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
                self.nextImage()
        else:
            # no file list found, load an image
            self.loadImage()

    def prevImage(self):
        if self.fileList:
            for _ in itertools.repeat(None, len(self.fileList) - 1):
                try:
                    self.filename = next(self.dirIterator)
                    
                except:
                    # the iterator has finished, restart it
                    self.dirIterator = iter(self.fileList)
                    self.nextImage()

            self.setWindowTitle(os.path.basename(self.filename))
            self.pixmap = QtGui.QPixmap(self.filename).scaled(self.label.size(), 
                QtCore.Qt.KeepAspectRatio)
            self.clearScene()
        else:
            self.loadImage()


    def clearScene(self):
        # print('clearing')
        self.label.scene.clear()
        self.label.BBList = []
        self.label.graphicsPixmapItem = QtWidgets.QGraphicsPixmapItem(QtGui.QPixmap(self.pixmap))
        self.label.scene.addItem(self.label.graphicsPixmapItem)

   

    def cropImage(self):
        if len(self.label.BBList) != 2:
            return
        savePath = os.path.join(os.getcwd(), 'Redo_Registration')
        # if not os.path.exists(savePath):
        #     os.mkdir(savePath)
        baseName = os.path.basename(self.filename)
        saveName = ''
        sInd = baseName.find('_s')
        if baseName[sInd - 3] == 'F':
            saveName = '1'
        else:
            saveName = '2'
        saveName = saveName + baseName[sInd - 1] + baseName[sInd + 2] + baseName[(sInd + 4):baseName.find('_error')] + '.txt'


        f = open(saveName, 'w')
        if self.label.BBList[0].top < self.label.BBList[1].top:
            f.write(str(self.label.BBList[0].top) + ' ' + str(self.label.BBList[0].bot) + ' ' + str(self.label.BBList[0].left) + ' ' + str(self.label.BBList[0].right)+ '\n')
            f.write(str(self.label.BBList[1].top) + ' ' + str(self.label.BBList[1].bot) + ' ' + str(self.label.BBList[1].left) + ' ' + str(self.label.BBList[1].right))
        else:

            f.write(str(self.label.BBList[1].top) + ' ' + str(self.label.BBList[1].bot) + ' ' + str(self.label.BBList[1].left) + ' ' + str(self.label.BBList[1].right)+ '\n')
            f.write(str(self.label.BBList[0].top) + ' ' + str(self.label.BBList[0].bot) + ' ' + str(self.label.BBList[0].left) + ' ' + str(self.label.BBList[0].right))
        
        f.close()





        



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    imageLoader = ImageLoader()
    imageLoader.show()
    sys.exit(app.exec_())