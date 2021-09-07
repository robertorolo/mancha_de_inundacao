#!/usr/bin/env python

"""gui.py: GUI do programa mancha_de_inundacao."""

__author__ = "Roberto Mentzingen Rolo"

#importando pacotes 
from tkinter import Tk, ttk, Label, Button, Entry, StringVar, Checkbutton, IntVar, INSERT
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
from tkinter.scrolledtext import ScrolledText
import geopandas
import rasterio
from rasterio.plot import show
import rasterio.mask
import matplotlib.pyplot as plt
import fiona
from mancha_de_inundacao import *
import ctypes
import logging
 
ctypes.windll.shcore.SetProcessDpiAwareness(1) #texto nitido em monitores de alta resolucao

fiona.supported_drivers['KML'] = 'rw'

class WidgetLogger(logging.Handler):
    def __init__(self, widget):
        logging.Handler.__init__(self)
        self.widget = widget

    def emit(self, record):
        # Append message (record) to the widget
        self.widget.insert(INSERT, record + '\n')
        self.widget.see("end")

def calcular_crio():
    longitude =  float(entry_x.get().replace(',','.'))
    latitude =  float(entry_y.get().replace(',','.'))
    global ponto_informado
    ponto_informado = Point((longitude, latitude))

    global h
    h = entry_h.get()
    h = float(h.replace(',','.'))

    global v
    v = entry_v.get()
    v = float(v.replace(',','.'))

    global criov
    criov = crio(v)

    wl1.emit('Comprimento do rio: {} km'.format(round(criov,2)))

    tabControl.tab(1,state="normal")

def carregar_tracado():
    wl2.emit('Lendo traçado do rio...')
    tracado_arquivo = askopenfilename()
    wl2.emit(tracado_arquivo)
    global tracado
    tracado = geopandas.read_file(tracado_arquivo, driver='KML')
    wl2.emit('Traçado do rio lido com sucesso!')
    btn_tracado["text"] = "Traçado do rio carregado"

def carregar_srtm():
    wl2.emit('Lendo SRTM...')
    srtm_arquivo = askopenfilename()
    wl2.emit(srtm_arquivo)
    global srtm
    srtm = rasterio.open(srtm_arquivo)
    wl2.emit('SRTM lido com sucesso!')
    btn_srtm["text"] = "SRTM carregado"

def calcular_perpendiculares():
    n = entry_suavi.get()
    n = int(n.replace(',','.'))

    comp = entry_comp.get()
    comp = float(comp.replace(',','.'))

    global tracado_simplificado
    tracado_simplificado = simplificar_tracado(tracado, n)

    global s
    global ds
    s, ds = secoes_perpendiculares(tracado_simplificado, n=21, comprimento=comp)
    s.crs = 'EPSG:31982'
    st = s.to_crs(epsg=4326)
    
    intersec_check = int(c_intersec_var.get())
    if intersec_check == 1:
        maxiter = int(entry_maxiter.get())
        maxtime = int(entry_maxtime.get())
        st = rotate_secs(st, wl=wl2, maxiter=maxiter, maxtime=maxtime, drange=[-10,10])

    wl2.emit('Feche a janela do mapa para contnuar.')
    fig, ax = plt.subplots(figsize=(8,8))
    
    show(srtm, ax=ax)
    st.plot(ax=ax, color='red')
    ax.scatter(ponto_informado.x, ponto_informado.y, color='red', label='Barragem')

    plt.legend()
    plt.show()

    if intersec_check == 1:
        st = st.to_crs(epsg=31982)
        s = st

    tabControl.tab(2,state="normal")
 
def calcular():
    wl3.emit('Cálculo hídrico iniciado...')

    qmax_barr = qmax_barragem(h, v)
    cotas(ponto_informado, srtm, h)

    c, dp, xs, ys = cotas_secoes(s, srtm)
    ct = [i[40] for i in c]

    alturas, qs = altura_de_agua_secoes(ds, dp, c, qmax_barr, v, h)

    x_all = []
    y_all = []
    h_all = []
    for idx, vv in enumerate(alturas):
        for idx1 in range(len(xs[idx])):
            h_all.append(vv)
            x_all.append(xs[idx][idx1])
            y_all.append(ys[idx][idx1])

    flname= 'srtm_cortado' #salva o srtm cortado. pode ser necessario apagar manualmente um arquivo salvo anteriormente

    clip_raster(s, srtm, flname)
    clipado = rasterio.open(flname)

    global xcoords
    global ycoords
    xcoords, ycoords, z = get_coordinates(clipado)
    #mascara
    chull = convex_hull(s)
    mascara = check_if_is_inside(chull, xcoords, ycoords)
    xcoords, ycoords, z = xcoords[mascara], ycoords[mascara], z[mascara]

    v_int = rbf_interpolation(x_all, y_all, h_all, xcoords, ycoords)
    global mancha
    mancha = np.where(v_int > z, 1, 0)

    wl3.emit('Feche a janela do mapa para continuar.')
    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(xcoords, ycoords, c=mancha, s=2)
    s.crs = 'EPSG:4326'
    st = s.to_crs(epsg=31982)
    s.plot(ax=ax, color='red', aspect=1)
    plt.show()

    wl3.emit('Finalizado!')
    entry_alpha['state'] = 'normal'
    entry_buffer['state'] = 'normal'
    btn_calculartab3['state'] = 'normal'
    
    for idx in range(len(qs)):
        wl3.emit('Seção {}: Vazão: {} - Altura da água: {} - Distância da barragem {}'.format(idx, round(qs[idx],2), round((alturas[idx]-ct[idx]),2), round(ds[idx],2)))

def importar_secoes():
    global s
    s = geopandas.read_file(shape_flname)
    btn_import["text"] = "Traçado importado"

def exportar_secoes():
    global shape_flname
    shape_flname = asksaveasfilename(defaultextension=".shp")
    
    exportar_geopandas(s, nome_do_arquivo=shape_flname)
    btn_export["text"] = "Traçado exportado"

def enable_maxiter():
    entry_maxiter['state'] = 'normal'
    entry_maxtime['state'] = 'normal'

def calcular_shape():

    alpha = entry_alpha.get()
    alpha = float(alpha.replace(',','.'))
    
    buffer = entry_buffer.get()
    buffer = float(buffer.replace(',','.'))

    global alpha_shape
    alpha_shape = mancha_pts_to_shape(xcoords, ycoords, mancha, alpha, buffer)
    
    wl3.emit('Feche a janela do mapa para continuar.')
    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(xcoords, ycoords, c=mancha, s=2)
    s.crs = 'EPSG:4326'
    st = s.to_crs(epsg=31982)
    s.plot(ax=ax, color='red', aspect=1)
    alpha_shape.plot(ax=ax, alpha=0.9, color='black')
    plt.show()
    
    btn_salvartab3['state'] = 'normal'

def salvartab3():
    kml_flname = asksaveasfilename()
    points_to_kml(xcoords, ycoords, mancha,kml_flname+str('_pts.kml'))
    alpha_shape.to_file(kml_flname+str('_shape.kml'), driver='KML')
    wl3.emit('Arquivos salvos!')

#GUI
root = Tk()
root.title("Mancha de inundação")
root.resizable(False, False)

#tab settings
tabControl = ttk.Notebook(root)

tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)
tab3 = ttk.Frame(tabControl)

tabControl.add(tab1, text='Passo 1')
tabControl.add(tab2, text='Passo 2', state="disabled")
tabControl.add(tab3, text='Passo 3',  state="disabled")

tabControl.pack(expand=1, fill="both")

#Tab1
#lat long
label_y = Label(tab1, text="Latitude:")
label_y.grid(row=1, column=0, sticky='W', padx=10, pady=10)

label_x = Label(tab1, text="Longitude:")
label_x.grid(row=2, column=0, sticky='W', padx=10, pady=10)

entry_y = Entry(tab1, width=8)
entry_y.grid(row=1, column=1, sticky='E', padx=10, pady=10)

entry_x = Entry(tab1, width=8)
entry_x.grid(row=2, column=1, sticky='E', padx=10, pady=10)

#h v
label_h = Label(tab1, text="Altura da barragem (m):")
label_h.grid(row=3, column=0, sticky='W', padx=10, pady=10)

label_v = Label(tab1, text="Volume armazenado (m³):")
label_v.grid(row=4, column=0, sticky='W', padx=10, pady=10)

entry_h = Entry(tab1, width=8)
entry_h.grid(row=3, column=1, sticky='E', padx=10, pady=10)

entry_v = Entry(tab1, width=8)
entry_v.grid(row=4, column=1, sticky='E', padx=10, pady=10)

#calcular
calcular_crio = Button(tab1, text="Calcular", command=calcular_crio)
calcular_crio.grid(row=5, column=1, sticky='E', padx=10, pady=10)

st1 = ScrolledText(tab1, height=10)
st1.grid(row=6, column=0, columnspan=2, sticky='E')
wl1 = WidgetLogger(st1)

#Tab2
btn_tracado = Button(tab2, text="Carregar traçado do rio", command=carregar_tracado)
btn_tracado.grid(row=0, column=0, sticky='W', padx=10, pady=10)

btn_srtm = Button(tab2, text="Carregar SRTM", command=carregar_srtm)
btn_srtm.grid(row=1, column=0, sticky='W', padx=10, pady=10)

label_suavi = Label(tab2, text="Número de retas para simplificação:")
label_suavi.grid(row=2, column=0, sticky='W', padx=10, pady=10)
entry_suavi = Entry(tab2, width=8)
entry_suavi.insert(0, "8")
entry_suavi.grid(row=2, column=1, sticky='E', padx=10, pady=10)

label_comp = Label(tab2, text="Comprimento das seções (m):")
label_comp.grid(row=3, column=0, sticky='W', padx=10, pady=10)
entry_comp = Entry(tab2, width=8)
entry_comp.insert(0, "4000")
entry_comp.grid(row=3, column=1, sticky='E', padx=10, pady=10)

#intersections
c_intersec_var = IntVar()
c_intersec = Checkbutton(tab2, text='Desinterceptar seções',variable=c_intersec_var, onvalue=1, offvalue=0, command=enable_maxiter)
c_intersec.grid(row=4, column=0, sticky='W', padx=10, pady=10)

label_maxiter = Label(tab2, text="Número máximo de iterações:")
label_maxiter.grid(row=5, column=0, sticky='W', padx=10, pady=10)
entry_maxiter = Entry(tab2, width=8)
entry_maxiter.insert(0, "1000")
entry_maxiter['state'] = 'disabled'
entry_maxiter.grid(row=5, column=1, sticky='E', padx=10, pady=10)

label_maxtime = Label(tab2, text="Tempo máximo (min):")
label_maxtime.grid(row=6, column=0, sticky='W', padx=10, pady=10)
entry_maxtime = Entry(tab2, width=8)
entry_maxtime.insert(0, "5")
entry_maxtime['state'] = 'disabled'
entry_maxtime.grid(row=6, column=1, sticky='E', padx=10, pady=10)

btn_export = Button(tab2, text="Exportar seções", command=exportar_secoes)
btn_export.grid(row=7, column=0, sticky='W', padx=10, pady=10)

btn_import = Button(tab2, text="Importar seções", command=importar_secoes)
btn_import.grid(row=7, column=1, sticky='E', padx=10, pady=10)

btn_calculartab2 = Button(tab2, text="Calcular", command=calcular_perpendiculares)
btn_calculartab2.grid(row=8, column=1, sticky='E', padx=10, pady=10)

st2 = ScrolledText(tab2, height=10)
st2.grid(row=9, column=0, columnspan=2, sticky='E')
wl2 = WidgetLogger(st2)

#tab3
label_alpha = Label(tab3, text="Alpha:")
label_alpha.grid(row=0, column=0, sticky='W', padx=10, pady=10)
entry_alpha = Entry(tab3, width=8)
entry_alpha.insert(0, "0.01")
entry_alpha['state'] = 'disabled'
entry_alpha.grid(row=0, column=1, sticky='E', padx=10, pady=10)

label_buffer = Label(tab3, text="Buffer:")
label_buffer.grid(row=1, column=0, sticky='W', padx=10, pady=10)
entry_buffer = Entry(tab3, width=8)
entry_buffer.insert(0, "30")
entry_buffer['state'] = 'disabled'
entry_buffer.grid(row=1, column=1, sticky='E', padx=10, pady=10)

btn_calculartab3 = Button(tab3, text="Calcular", command=calcular)
btn_calculartab3.grid(row=2, column=0, sticky='W', padx=10, pady=10)

btn_calculartab3 = Button(tab3, text="(Re)calcular shape", command=calcular_shape)
btn_calculartab3['state'] = 'disabled'
btn_calculartab3.grid(row=2, column=1, sticky='E', padx=10, pady=10)

btn_salvartab3 = Button(tab3, text="Salvar", command=salvartab3)
btn_salvartab3['state'] = 'disabled'
btn_salvartab3.grid(row=3, column=0, sticky='W', padx=10, pady=10)

st3 = ScrolledText(tab3, height=10)
st3.grid(row=4, column=0, columnspan=2, sticky='E')
wl3 = WidgetLogger(st3)

root.mainloop()