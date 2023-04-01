from PyQt5 import QtCore, QtGui, QtWidgets
from mancha_de_inundacao import *

import fiona
fiona.supported_drivers['KML'] = 'rw'

import matplotlib.pyplot as plt
from rasterio.plot import show
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from pyvistaqt import QtInteractor
import pyvista as pv

import pandas as pd
import random

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        
        self.fig, self.ax = plt.subplots(figsize=(10,10))
        self.canvas = FigureCanvas(self.fig)
        
        self.plotter =QtInteractor()
    
        Dialog.setObjectName("Dialog")
        Dialog.resize(900, 900)
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.gridLayout.setObjectName("gridLayout")
        self.tabWidget = QtWidgets.QTabWidget(Dialog)
        self.tabWidget.setEnabled(True)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_3 = QtWidgets.QLabel(self.tab)
        self.label_3.setObjectName("label_3")
        self.gridLayout_3.addWidget(self.label_3, 2, 0, 1, 1)
        self.hval = QtWidgets.QDoubleSpinBox(self.tab)
        self.hval.setMaximum(100000000000.0)
        self.hval.setObjectName("hval")
        self.gridLayout_3.addWidget(self.hval, 2, 2, 1, 1)
        self.label = QtWidgets.QLabel(self.tab)
        self.label.setObjectName("label")
        self.gridLayout_3.addWidget(self.label, 0, 0, 1, 1)
        self.longval = QtWidgets.QDoubleSpinBox(self.tab)
        self.longval.setMaximum(1000000.0)
        self.longval.setMinimum(-1000000.0)
        self.longval.setDecimals(4)
        self.longval.setObjectName("longval")
        self.gridLayout_3.addWidget(self.longval, 1, 2, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.tab)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 3, 0, 1, 1)
        self.latval = QtWidgets.QDoubleSpinBox(self.tab)
        self.latval.setMaximum(1000000.0)
        self.latval.setMinimum(-1000000.0)
        self.latval.setDecimals(4)
        self.latval.setObjectName("latval")
        self.gridLayout_3.addWidget(self.latval, 0, 2, 1, 1)
        self.tab1calc = QtWidgets.QPushButton(self.tab)
        self.tab1calc.setObjectName("tab1calc")
        self.gridLayout_3.addWidget(self.tab1calc, 4, 2, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.tab)
        self.label_2.setObjectName("label_2")
        self.gridLayout_3.addWidget(self.label_2, 1, 0, 1, 1)
        self.vval = QtWidgets.QDoubleSpinBox(self.tab)
        self.vval.setMaximum(100000000000.0)
        self.vval.setObjectName("vval")
        self.gridLayout_3.addWidget(self.vval, 3, 2, 1, 1)
        self.crioval = QtWidgets.QLabel(self.tab)
        self.crioval.setObjectName("crioval")
        self.gridLayout_3.addWidget(self.crioval, 5, 2, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.tab)
        self.label_8.setObjectName("label_8")
        self.gridLayout_3.addWidget(self.label_8, 5, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem, 6, 0, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.tab_2)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.nretasval = QtWidgets.QSpinBox(self.tab_2)
        self.nretasval.setEnabled(True)
        self.nretasval.setProperty("value", 8)
        self.nretasval.setObjectName("nretasval")
        self.gridLayout_4.addWidget(self.nretasval, 2, 1, 1, 1)
        self.matplotlibwidget = self.canvas
        self.matplotlibwidget.setObjectName("matplotlibwidget")
        self.gridLayout_4.addWidget(self.matplotlibwidget, 8, 0, 1, 2)
        self.mtoolbar = NavigationToolbar(self.canvas, self.matplotlibwidget)
        self.tab2calc = QtWidgets.QPushButton(self.tab_2)
        self.tab2calc.setEnabled(True)
        self.tab2calc.setObjectName("tab2calc")
        self.gridLayout_4.addWidget(self.tab2calc, 7, 1, 1, 1)
        self.carregartracado = QtWidgets.QPushButton(self.tab_2)
        self.carregartracado.setEnabled(True)
        self.carregartracado.setObjectName("carregartracado")
        self.gridLayout_4.addWidget(self.carregartracado, 1, 1, 1, 1)
        self.carregarsrtm = QtWidgets.QPushButton(self.tab_2)
        self.carregarsrtm.setEnabled(True)
        self.carregarsrtm.setObjectName("carregarsrtm")
        self.gridLayout_4.addWidget(self.carregarsrtm, 0, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.tab_2)
        self.label_5.setEnabled(True)
        self.label_5.setObjectName("label_5")
        self.gridLayout_4.addWidget(self.label_5, 2, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.tab_2)
        self.label_6.setEnabled(True)
        self.label_6.setObjectName("label_6")
        self.gridLayout_4.addWidget(self.label_6, 3, 0, 1, 1)
        self.exportsec = QtWidgets.QPushButton(self.tab_2)
        self.exportsec.setEnabled(True)
        self.exportsec.setObjectName("exportsec")
        self.gridLayout_4.addWidget(self.exportsec, 5, 1, 1, 1)
        self.comsecval = QtWidgets.QSpinBox(self.tab_2)
        self.comsecval.setEnabled(True)
        self.comsecval.setMaximum(1000000)
        self.comsecval.setProperty("value", 4000)
        self.comsecval.setObjectName("comsecval")
        self.gridLayout_4.addWidget(self.comsecval, 3, 1, 1, 1)
        self.importsec = QtWidgets.QPushButton(self.tab_2)
        self.importsec.setEnabled(True)
        self.importsec.setObjectName("importsec")
        self.gridLayout_4.addWidget(self.importsec, 6, 1, 1, 1)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.tab_3)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.saveshape = QtWidgets.QPushButton(self.tab_3)
        self.saveshape.setObjectName("saveshape")
        self.gridLayout_2.addWidget(self.saveshape, 1, 0, 1, 1)
        self.tab3calc = QtWidgets.QPushButton(self.tab_3)
        self.tab3calc.setObjectName("tab3calc")
        self.gridLayout_2.addWidget(self.tab3calc, 0, 0, 1, 1)
        self.savereport = QtWidgets.QPushButton(self.tab_3)
        self.savereport.setObjectName("savereport")
        self.gridLayout_2.addWidget(self.savereport, 2, 0, 1, 1)
        self.pyvistawidget = self.plotter
        self.pyvistawidget.setObjectName("pyvistawidget")
        self.gridLayout_2.addWidget(self.pyvistawidget, 3, 0, 1, 1)
        self.tabWidget.addTab(self.tab_3, "")
        self.gridLayout.addWidget(self.tabWidget, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

        Dialog.setTabOrder(self.latval, self.longval)
        Dialog.setTabOrder(self.longval, self.hval)
        Dialog.setTabOrder(self.hval, self.vval)
        Dialog.setTabOrder(self.vval, self.tab1calc)

        Dialog.setTabOrder(self.tab1calc, self.carregarsrtm)
        Dialog.setTabOrder(self.carregarsrtm, self.carregartracado)
        Dialog.setTabOrder(self.carregartracado, self.nretasval)
        Dialog.setTabOrder(self.nretasval, self.comsecval)
        Dialog.setTabOrder(self.comsecval, self.exportsec)
        Dialog.setTabOrder(self.exportsec, self.importsec)
        Dialog.setTabOrder(self.importsec, self.tab2calc)

        Dialog.setTabOrder(self.tab2calc, self.tab3calc)
        Dialog.setTabOrder(self.tab3calc, self.saveshape)
        Dialog.setTabOrder(self.saveshape, self.savereport)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Mancha de inundação"))
        self.label_3.setText(_translate("Dialog", "Altura da barragem (m)"))
        self.label.setText(_translate("Dialog", "Latitude"))
        self.label_4.setText(_translate("Dialog", "Volume armazenado (m³)"))
        self.tab1calc.setText(_translate("Dialog", "Calcular"))
        self.label_2.setText(_translate("Dialog", "Longitude"))
        self.crioval.setText(_translate("Dialog", "0"))
        self.label_8.setText(_translate("Dialog", "Comprimento do rio (km)"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("Dialog", "Passo 1"))
        self.tab2calc.setText(_translate("Dialog", "Calcular"))
        self.carregartracado.setText(_translate("Dialog", "Carregar traçado"))
        self.carregarsrtm.setText(_translate("Dialog", "Carregar SRTM"))
        self.label_5.setText(_translate("Dialog", "Número de retas para simplificação"))
        self.label_6.setText(_translate("Dialog", "Comprimento da seção (m)"))
        self.exportsec.setText(_translate("Dialog", "Exportar seções"))
        self.importsec.setText(_translate("Dialog", "Importar seções"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("Dialog", "Passo 2"))
        self.saveshape.setText(_translate("Dialog", "Salvar shape file"))
        self.tab3calc.setText(_translate("Dialog", "Calcular"))
        self.savereport.setText(_translate("Dialog", "Salvar relatório"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("Dialog", "Passo 3"))
        
        self.tab1calc.clicked.connect(self.calcrio)
        
        self.carregarsrtm.clicked.connect(self.carregar_srtm)
        self.carregartracado.clicked.connect(self.carregar_tracado)
        self.tab2calc.clicked.connect(self.calcular_perpendiculares)
        self.exportsec.clicked.connect(self.exportar_secoes)
        self.importsec.clicked.connect(self.importar_secoes)
        
        self.saveshape.clicked.connect(self.save_shp)
        self.savereport.clicked.connect(self.save_report)
        self.tab3calc.clicked.connect(self.calcular)
        
    def calcrio(self):
        print('Calculando comprimento do rio...')
        longitude = self.longval.value()
        latitude =  self.latval.value()
  
        self.ponto_informado = Point((longitude, latitude))

        self.h = self.hval.value()
        self.v = self.vval.value()

        self.criov = crio(self.v)
        self.criov = round(self.criov, 2)

        self.crioval.setText(str(self.criov))
        
    def carregar_srtm(self):
        file_path = QtWidgets.QFileDialog.getOpenFileName(None, "Selecione o SRTM", "", "tif files (*.tif)")
        srtm_arquivo = file_path[0] 
        self.srtm = rasterio.open(srtm_arquivo)
        print('SRTM carregado!')
        
    def carregar_tracado(self):
        file_path = QtWidgets.QFileDialog.getOpenFileName(None, "Selecione o traçado", "", "kml files (*.kml)")
        tracado_arquivo = file_path[0] 
        self.tracado = geopandas.read_file(tracado_arquivo, driver='KML')
        print('Traçado carregado!')

    def calcular_perpendiculares(self):
        print('Calculando perpendiculares...')
        self.ax.clear()
        n = self.nretasval.value()
        comp = self.comsecval.value()

        self.tracado_simplificado = simplificar_tracado(self.tracado, n)

        self.s, self.ds = secoes_perpendiculares(self.tracado_simplificado, n=21, comprimento=comp)
        self.s.crs = f'EPSG:{datum}'
        st = self.s.to_crs(epsg=4326)
        
        show(self.srtm, ax=self.ax)
        st.plot(ax=self.ax, color='red')
        self.ax.scatter(self.ponto_informado.x, self.ponto_informado.y, color='red', label='Barragem')
        
    def exportar_secoes(self):
        self.shape_flname = QtWidgets.QFileDialog.getSaveFileName(None, "Selecione onde o arquivo será salvo", "", "shp files (*.shp)")
        self.shape_flname = self.shape_flname[0]
        exportar_geopandas(self.s, nome_do_arquivo=self.shape_flname)
        
    def importar_secoes(self):
        self.s = geopandas.read_file(self.shape_flname)
        
        self.ax.clear()
        show(self.srtm, ax=self.ax)
        self.s.plot(ax=self.ax, color='red')
        self.ax.scatter(self.ponto_informado.x, self.ponto_informado.y, color='red', label='Barragem')
        
        self.s = self.s.to_crs(epsg=datum)
        
    def save_shp(self):
        self.shape_flname_mancha = QtWidgets.QFileDialog.getSaveFileName(None, "Selecione onde o arquivo será salvo", "", "shp files (*.shp)")
        flname = self.shape_flname_mancha[0]
        surfaces_to_kml(self.surf_surface, self.surf_water, flname)
        print('Shape file salvo!')

    def save_report(self):
        nomes = ['seção {}'.format(i) for i in range(21)]
        alturas_a = [self.alturas[idx]-self.ct[idx] for idx in range(21)]
        data_array = np.array(
            [nomes,
            self.ds,
            self.qs,
            alturas_a]
        ).T
        df = pd.DataFrame(columns=['Seções', 'Distância', 'Vazão', 'Altura de água'], data=data_array)
        self.rela_flname = QtWidgets.QFileDialog.getSaveFileName(None, "Selecione onde o arquivo será salvo", "", "CSV files (*.csv)")
        flname = self.rela_flname[0]
        df.to_csv(flname, index=False)
        print('Relatório salvo!')

    def calcular(self):
        print('Calculo hídrico iniciado...')
        fc = 1

        qmax_barr = qmax_barragem(self.h, self.v)
        cotas(self.ponto_informado, self.srtm, self.h)

        c, dp, xs, ys = cotas_secoes(self.s, self.srtm)
        self.ct = [i[40] for i in c]

        self.alturas, self.qs = altura_de_agua_secoes(self.ds, dp, c, qmax_barr, self.v, self.h, fc)

        x_all = []
        y_all = []
        h_all = []
        for idx, vv in enumerate(self.alturas):
            for idx1 in range(len(xs[idx])):
                h_all.append(vv)
                x_all.append(xs[idx][idx1])
                y_all.append(ys[idx][idx1])

        int1 = random.randrange(0,9,1)
        int2 = random.randrange(0,9,1)
        flname= 'srtm_cortado_'+str(int1)+str(int2)

        clip_raster(self.s, self.srtm, flname)
        clipado = rasterio.open(flname)

        xcoords, ycoords, z = get_coordinates(clipado)

        chull = convex_hull(self.s)
        mascara = check_if_is_inside(chull, xcoords, ycoords)
        xcoords, ycoords, z = xcoords[mascara], ycoords[mascara], z[mascara]

        pts = [[i, j, k] for i, j, k in zip(x_all, y_all, h_all)]
        pts = np.array(pts)
        pts_surf = [[i,j,k] for i, j, k in zip(xcoords, ycoords, z)]
        pts_surf = np.array(pts_surf)
       
        cloud1 = pv.PolyData(pts)
        cloud2 = pv.PolyData(pts_surf)

        self.surf_water = cloud1.delaunay_2d()
        self.surf_surface = cloud2.delaunay_2d()
        self.plotter.add_mesh(self.surf_water, name='water surface', color='blue')
        self.plotter.add_mesh(self.surf_surface, name='terrain surface', cmap='viridis', scalars=z)
        self.plotter.view_isometric()

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
