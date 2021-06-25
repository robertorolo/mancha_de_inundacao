#importando pacotes 
from tkinter import Tk, Label, Button, Entry
from tkinter.filedialog import askopenfilename, asksaveasfile
import geopandas
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
import fiona
from shapely.geometry import Point
import numpy as np

fiona.supported_drivers['KML'] = 'rw'

#definindo funcoes
def crio(volume):
    volume = volume * 10e-6
    # Calcula o comprimento do rio a ser modelado (sobre o rio suavizado)
    crio = 0.0000000887*float(volume)**3 - 0.00026*float(volume)**2 + 0.265*float(volume) + 6.74
    if crio < 5.0:
        crio = 5.0
    if crio > 100.0:
        crio = 100.0
    return crio

def qmax_barragem(altura, volume):
    mmc = 0.0039 * volume ** 0.8122
    froe = 0.607 * (volume ** 0.295 * altura ** 1.24)
    return max(mmc, froe)

def qmax_secao(x, q_max_barr, volume):
    if volume * 10e-6 > 6.2:
        return q_max_barr * 10 ** (-0.01243*x)
    else:
        a = 0.002 * np.log(volume) + 0.9626
        b = -0.20047 * (volume + 25000) ** -0.5979
        return q_max_barr * a * np.exp(b*x)

def suavizar_tracado():
    print(tracado)
    tracado_simplificado = tracado['geometry'].simplify(5)

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

def calcular1():
    longitude =  float(entry_x.get().replace(',','.'))
    latitude =  float(entry_y.get().replace(',','.'))
    global ponto_informado
    ponto_informado = Point((longitude, latitude))
    h =  float(entry_h.get().replace(',','.'))
    v =  float(entry_v.get().replace(',','.'))
    global criov
    criov = crio(v)
    print('Comprimento do rio {} km'.format(criov))
    #global cota_tal
    #cota_tal = list(srtm.sample([(ponto_informado.x, ponto_informado.y)]))
    #cota_tal = cota_tal[0][0]
    #global cota_cor 
    #cota_cor = cota_tal + h
    #print('cota do talvegue: {}; cota do coroamento {}'.format(cota_tal, cota_cor))
    global q_barr 
    q_barr = qmax_barragem(h, v)
    print('vazão máxima na seção da barragem: {} m³/s'.format(q_barr))

def calcular2():
    pass

def plotar2():
    fig, ax = plt.subplots(figsize=(8,8))
    show(srtm, ax=ax)

    ax.scatter(ponto_informado.x, ponto_informado.y, color='red', label='Barragem')

    tracado.plot(ax=ax, color='red')
    
    plt.legend()
    plt.show()
    
#GUI
root = Tk()
root.title("Mancha de inundação")
root.resizable(False, False)

#label passos
label_p1 = Label(root, text="Passo 1")
label_p1.grid(row=0, column=0, columnspan=2, padx=10, pady=10)

label_p2 = Label(root, text="Passo 2")
label_p2.grid(row=0, column=2, columnspan=2, padx=10, pady=10)

#lat long
label_y = Label(root, text="Latitude:")
label_y.grid(row=1, column=0, sticky='W', padx=10, pady=10)

label_x = Label(root, text="Longitude:")
label_x.grid(row=2, column=0, sticky='W', padx=10, pady=10)

entry_y = Entry(root, width=8)
entry_y.grid(row=1, column=1, sticky='E', padx=10, pady=10)

entry_x = Entry(root, width=8)
entry_x.grid(row=2, column=1, sticky='E', padx=10, pady=10)

#h v
label_h = Label(root, text="Altura da barragem (m):")
label_h.grid(row=3, column=0, sticky='W', padx=10, pady=10)

label_v = Label(root, text="Volume armazenado (m³):")
label_v.grid(row=4, column=0, sticky='W', padx=10, pady=10)

entry_h = Entry(root, width=8)
entry_h.grid(row=3, column=1, sticky='E', padx=10, pady=10)

entry_v = Entry(root, width=8)
entry_v.grid(row=4, column=1, sticky='E', padx=10, pady=10)

#btns
btn_tracado = Button(root, text="Carregar traçado do rio", command=carregar_tracado)
btn_tracado.grid(row=1, column=2, sticky='W', padx=10, pady=10)

btn_srtm = Button(root, text="Carregar SRTM", command=carregar_srtm)
btn_srtm.grid(row=2, column=2, sticky='W', padx=10, pady=10)

btn_calcular1 = Button(root, text="Calcular", command=calcular1)
btn_calcular1.grid(row=7, column=1, sticky='E', padx=10, pady=10)

btn_calcular2 = Button(root, text="Calcular", command=calcular2)
btn_calcular2.grid(row=3, column=2, sticky='E', padx=10, pady=10)

btn_plotar2 = Button(root, text="Plotar", command=plotar2)
btn_plotar2.grid(row=4, column=2, sticky='E', padx=10, pady=10)

root.mainloop()