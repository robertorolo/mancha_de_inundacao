import numpy as np
from shapely.geometry import Point, LineString
from pyproj import Proj, transform

def crio(volume):
    volume = volume * 10e-6 # tranformacao para hm
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

    return tracado

def exportar_geopandas(tracado, nome_do_arquivo='tracado.shp'):
    tracado.crs = 'EPSG:31982'
    tracado = tracado.to_crs(epsg=4674)
    tracado.to_file(nome_do_arquivo)

def cotas_secoes(tracado, srtm):
    inProj = Proj('epsg:31982')
    outProj = Proj('epsg:4674')
    cotas = []
    for linha in tracado.iloc[2:]['geometry']:
        p, d = pontos_tracado(linha, n=81)
        p = [i.coords[:][0] for i in p]
        pt = [(transform(inProj,outProj,i[0],i[1])) for i in p]
        pt = [(i[1], i[0]) for i in pt]
        cota = list(srtm.sample(pt))
        cota = [k[0] for k in cota]
        cotas.append(cota)
    return cotas, d
