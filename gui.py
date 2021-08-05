#!/usr/bin/env python

"""gui.py: GUI do programa mancha_de_inundacao."""

__author__ = "Roberto Mentzingen Rolo"

#importando pacotes 
from tkinter import Tk, ttk, Label, Button, Entry, StringVar
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
import geopandas
import rasterio
from rasterio.plot import show
import rasterio.mask
import matplotlib.pyplot as plt
import fiona
from mancha_de_inundacao import *
import ctypes
 
ctypes.windll.shcore.SetProcessDpiAwareness(1) #texto nitido em monitores de alta resolucao

fiona.supported_drivers['KML'] = 'rw'

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
    crio_var.set(str(round(criov,2)))

def carregar_tracado():
    print('Lendo traçado do rio...')
    tracado_arquivo = askopenfilename()
    print(tracado_arquivo)
    global tracado
    tracado = geopandas.read_file(tracado_arquivo, driver='KML')
    print('Traçado do rio lido com sucesso!\n')
    btn_tracado["text"] = "Traçado do rio carregado"

def carregar_srtm():
    print('Lendo SRTM...')
    srtm_arquivo = askopenfilename()
    print(srtm_arquivo)
    global srtm
    srtm = rasterio.open(srtm_arquivo)
    print('SRTM lido com sucesso!\n')
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

    print('Feche a janela do mapa para contnuar.')
    fig, ax = plt.subplots(figsize=(8,8))
    
    show(srtm, ax=ax)
    st.plot(ax=ax, color='red')
    ax.scatter(ponto_informado.x, ponto_informado.y, color='red', label='Barragem')

    plt.legend()
    plt.show()
 
def calcular():
    print('Cálculo hídrico iniciado...')

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

    xcoords, ycoords, z = get_coordinates(clipado)

    v_int = rbf_interpolation(x_all, y_all, h_all, xcoords, ycoords)
    mancha = np.where(v_int > z, 1, 0)
    
    print('Feche a janela do mapa para continuar.')
    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(xcoords, ycoords, c=mancha, s=2)
    s.crs = 'EPSG:4326'
    st = s.to_crs(epsg=31982)
    s.plot(ax=ax, color='red')
    plt.show()

    kml_flname = asksaveasfilename(defaultextension=".kml")
    points_to_kml(xcoords, ycoords, mancha, kml_flname)

    print('Finalizado!')
    
    for idx in range(len(qs)):
        print('Seção {}: Vazão: {} - Altura da água: {} - Distância da barragem {}'.format(idx, round(qs[idx],2), round((alturas[idx]-ct[idx]),2), round(ds[idx],2)))

def importar_secoes():
    global s
    s = geopandas.read_file(shape_flname)
    btn_import["text"] = "Traçado importado"

def exportar_secoes():
    global shape_flname
    shape_flname = asksaveasfilename(defaultextension=".shp")
    
    exportar_geopandas(s, nome_do_arquivo=shape_flname)
    btn_export["text"] = "Traçado exportado"

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
tabControl.add(tab2, text='Passo 2')
tabControl.add(tab3, text='Passo 3')

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

#crio
label_crio = Label(tab1, text="Comprimento do rio (km):")
label_crio.grid(row=5, column=0, sticky='W', padx=10, pady=10)

crio_var = StringVar()
label_crio_resultado = Label(tab1, textvariable=crio_var)
label_crio_resultado.grid(row=5, column=1, sticky='E', padx=10, pady=10)

#calcular
calcular_crio = Button(tab1, text="Calcular", command=calcular_crio)
calcular_crio.grid(row=6, column=1, sticky='E', padx=10, pady=10)

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

btn_export = Button(tab2, text="Exportar seções", command=exportar_secoes)
btn_export.grid(row=4, column=0, sticky='W', padx=10, pady=10)

btn_import = Button(tab2, text="Importar seções", command=importar_secoes)
btn_import.grid(row=4, column=1, sticky='W', padx=10, pady=10)

btn_calculartab2 = Button(tab2, text="Calcular", command=calcular_perpendiculares)
btn_calculartab2.grid(row=5, column=1, sticky='E', padx=10, pady=10)

btn_calculartab3 = Button(tab3, text="Calcular", command=calcular)
btn_calculartab3.grid(row=0, column=0, sticky='E', padx=10, pady=10)

root.mainloop()