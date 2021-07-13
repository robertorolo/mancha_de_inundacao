import numpy as np
import numpy.polynomial.polynomial as poly
from shapely.geometry import Point, LineString, box
import pyproj
from pyproj import Proj, Transformer, transform
import math
import rasterio.mask
import rasterio

def crio(volume):
    volume = volume * 10e-6 # tranformacao para hm
    # Calcula o comprimento do rio a ser modelado (sobre o rio suavizado)
    crio = 0.0000000887*float(volume)**3 - 0.00026*float(volume)**2 + 0.265*float(volume) + 6.74
    if crio < 5.0:
        crio = 5.0
    elif crio > 100.0:
        crio = 100.0
    else:
        crio = crio
    return crio

def qmax_barragem(altura, volume):
    mmc = 0.0039 * volume ** 0.8122
    froe = 0.607 * (volume ** 0.295 * altura ** 1.24)
    return max(mmc, froe)

def qmax_secao(x, q_max_barr, volume):
    volume = volume * 10e-6
    x = x * 10e-3
    if volume > 6.2:
        return q_max_barr * 10 ** (-0.01243*x)
    else:
        a = 0.002 * np.log(volume) + 0.9626
        b = -0.20047 * (volume + 25000) ** -0.5979
        return q_max_barr * a * np.exp(b*x)

def pontos_tracado(linha, n=21):
    distances = np.linspace(0, linha.length, n)
    points = [linha.interpolate(distance) for distance in distances]

    return points, distances

def cotas(ponto_informado, srtm, altura):
    cota_tal = list(srtm.sample([(ponto_informado.x, ponto_informado.y)]))
    cota_tal = cota_tal[0][0]

    cota_cor = cota_tal + altura

    return cota_tal, cota_cor

def simplificar_tracado(tracado, n):
    tracado = tracado.to_crs(epsg=31982)

    distances = np.linspace(0, tracado.length, n)
    points = [tracado.interpolate(distance) for distance in distances]

    x = [point.x[0] for point in points]
    y = [point.y[0] for point in points]

    ls = LineString([Point(i, j) for i, j in zip(x, y)])
    data = ['tracado do rio simplificado', '', ls]
    tracado.loc[len(tracado)] = data

    return tracado

def split_linha(tracado):
    line_split = []
    coords = tracado.iloc[1]['geometry'].coords[:]
    for i in range(len(coords)):
        if i < len(coords)-1:
            a = coords[i]
            b = coords[i+1]
            l = LineString([a, b])
            line_split.append(l)

    return line_split

def perpendicular(linha, ponto, comprimento=4000):
    rp1, rp2 = linha.coords[:][0], linha.coords[:][1]
    slope=(rp2[1]-rp1[1])/(rp2[0]-rp1[0])
    (rp2[1]-rp1[1])/(rp2[0]-rp1[0])

    C, D = [0,0], [0,0]
    B = ponto.coords[:][0]
    dy = np.sqrt((comprimento/2)**2/(slope**2+1))
    dx = -slope*dy
    C[0] = B[0] + dx
    C[1] = B[1] + dy
    D[0] = B[0] - dx
    D[1] = B[1] - dy
    perp = LineString([tuple(C), tuple(D)])

    return perp

def secoes_perpendiculares(tracado, n=21, comprimento=4000):
    p, d = pontos_tracado(tracado.iloc[1]['geometry'], n)
    l = split_linha(tracado)
    tol = 1e-8
    for i, point in enumerate(p):
        for line in l:
            if line.distance(point) < tol:
                perp = perpendicular(line, point)
                data = ['seçao {}'.format(i), d[i], perp]
                tracado.loc[len(tracado)] = data

    return tracado, d

def exportar_geopandas(tracado, nome_do_arquivo='tracado.shp'):
    tracado.crs = 'EPSG:31982'
    tracado = tracado.to_crs(epsg=4326)
    tracado.to_file(nome_do_arquivo)

def transformacao(x, y, d_to_m, new):#takes a loot of time to run
    if d_to_m:
        outProj = "epsg:31982"
        inProj = "epsg:4326"
    else:
        inProj = "epsg:31982"
        outProj = "epsg:4326"

    if new:
        transformer = Transformer.from_crs(inProj, outProj)
        np = transformer.transform(x, y)
    else:
        np = transform(inProj, outProj, x, y)
    return np

def cotas_secoes(tracado, srtm):
    cotas = []
    for linha in tracado.iloc[2:]['geometry']:
        p, d = pontos_tracado(linha, n=81)
        p = [i.coords[:][0] for i in p]
        x, y = [i[0] for i in p], [j[1] for j in p]
        pt = transformacao(x, y, d_to_m=False, new=True)
        pt = [(y, x) for x,y in zip(pt[0], pt[1])]
        cota = list(srtm.sample(pt))
        cota = [k[0] for k in cota]
        cotas.append(cota)
    return cotas, d

def raio_hidraulico(y, x):
    y, x = np.array(y), np.array(x)
    yt = 1 - y + max(y)
    hs = np.linspace(0, max(yt), 11)
    areas = []
    radius = []
    for h in hs[1:]:
        ytt = yt - (max(yt) - h)
        f = ytt > 0
        ytt, xt = ytt[f], x[f]
        area = np.trapz(y=ytt, x=xt)
        distances = []
        for i in range(len(ytt)):
            if i < len(ytt) - 1:
                d = math.hypot(xt[i+1]-xt[i], ytt[i+1]-ytt[i])
                distances.append(d)
        perimeter = np.sum(distances)
        areas.append(area)
        radius.append(area/perimeter)
    return areas, radius, hs

def manning(a, r, j, k=15):
    a, r = np.array(a), np.array(r)
    q = k*a*r**(2/3)*j**(1/2)
    return q

def polyfit(x, y, x_i):
    coefs = poly.polyfit(x, y, 3)
    ffit = poly.polyval(x_i, coefs)
    return ffit

def clip_raster(secs, srtm):
    secs.crs = 'EPSG:31982'
    secs = secs.to_crs(epsg=4326)
    minx, miny, maxx, maxy = min(secs.bounds['minx']), min(secs.bounds['miny']), max(secs.bounds['maxx']), max(secs.bounds['maxy'])
    bbox = box(minx, miny, maxx, maxy)
    out_image, out_transform = rasterio.mask.mask(srtm, [bbox], crop=True)
    out_meta = srtm.meta
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})
    with rasterio.open('srtm_cropado', 'w', **out_meta) as dest:
        dest.write(out_image)

def get_coordinates(clipado):
    w = clipado.width
    h = clipado.height
    ij = []
    for i in range(h):
        for j in range(w):
            a = clipado.xy(i, j, offset='center')
            ij.append(a)
    b = clipado.read(1).flatten()
    x, y = [i[0] for i in ij], [j[1] for j in ij]
    ij = transformacao(y, x, d_to_m=True, new=True)
    return ij, b
