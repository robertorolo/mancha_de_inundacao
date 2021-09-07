"""mancha_de_inundacao.py: modulo com as funcoes necessarias para o programa mancha_de_inundacao."""

__author__ = "Roberto Mentzingen Rolo"

#importando pacotes
import numpy as np
import numpy.polynomial.polynomial as poly
from shapely.geometry import Point, LineString, box, MultiLineString
from shapely import affinity
from pyproj import Transformer, transform
import math
import rasterio.mask
import rasterio
import pandas as pd
import geopandas
from time import time

def check_if_is_inside(chull, x, y):
    #retorna uma mascara pertence ao poligono
    mask = []
    pts = [Point(i, j) for i,j in zip(x, y)]
    for i in pts:
        if chull.contains(i):
            mask.append(True)
        else:
            mask.append(False)
    return np.array(mask)

def convex_hull(s):
    #area em torno das secoes
    a = geopandas.GeoSeries([i for i in s.iloc[2:].geometry])
    b = a.unary_union.convex_hull
    return b

def rotate_l(l1, drange):
    #gira um linha em seu centro um valor aleatorio em graus entre um range definido
    g = np.random.uniform(drange[0],drange[1])
    nl = affinity.rotate(l1, g, 'center')
    return nl

def check_if_intercepts(l1, l2):
    #verifica se duas secoes se interceptam
    return l1.intersects(l2)

def check_all(ndf):
    #verifica se uma secao se intercepta com qualquer outra
    intersect = False
    for i in range(ndf.shape[0]-1):
        for j in range(ndf.shape[0]-1):
            if i != j:
                l1 = ndf.iloc[i].geometry
                l2 = ndf.iloc[j].geometry
                if check_if_intercepts(l1, l2):
                    intersect = True
    return intersect

def rotate_secs(sec_df, wl, maxiter=1000, maxtime=5, drange=[-10,10]):
    #tenta desiterceptar as seções
    t1 = time()
    delta_t = 0
    ndf = sec_df.iloc[2:].copy(deep=False)
    ndf = ndf.reset_index(drop=True)

    while check_all(ndf) and delta_t <= maxtime*60:
        t2 = time()
        delta_t = t2-t1
    
        for idx in range(ndf.shape[0]-1):
            if idx == 0:
                l1 = ndf.iloc[idx].geometry
                l2 = ndf.iloc[idx+1].geometry
                itera = 0
                while check_if_intercepts(l1, l2) and itera <= maxiter:
                    l1 = rotate_l(l1, drange)
                    itera = itera + 1
                ndf.iloc[idx].geometry = l1

            elif idx == ndf.shape[0]-1:
                l1 = ndf.iloc[idx].geometry
                l2 = ndf.iloc[idx-1].geometry
                itera = 0
                while check_if_intercepts(l1, l2) and itera <= maxiter:
                    l1 = rotate_l(l1, drange)
                    itera = itera + 1
                ndf.geometry.iloc[idx] = l1

            else:
                l1 = ndf.iloc[idx].geometry
                l2 = ndf.iloc[idx-1].geometry
                l3 = ndf.iloc[idx+1].geometry
                itera = 0
                while (check_if_intercepts(l1, l2) or check_if_intercepts(l1, l3)) and itera <= maxiter:
                    l1 = rotate_l(l1, drange)
                    itera = itera + 1
                ndf.geometry.iloc[idx] = l1

    ndf = geopandas.GeoDataFrame(pd.concat([sec_df.iloc[:2], ndf], ignore_index=True))

    if delta_t >= maxtime:
        wl.emit('Nao foi possivel desinterceptar as secoes apos {} minutos. Voce pode tentar novamente.'.format(maxtime))
    else:
        wl.emit('Isto levou {} segundos.'.format(int(delta_t)))
    
    return ndf

def crio(volume):
    #calcula o comprimento do rio a ser modelado a partir da barragem
    volume = volume * 1E-06 # tranformacao para hm
    crio = 0.0000000887*float(volume)**3 - 0.00026*float(volume)**2 + 0.265*float(volume) + 6.74
    if crio < 5.0:
        crio = 5.0
    elif crio > 100.0:
        crio = 100.0
    else:
        crio = crio
    return crio

def qmax_barragem(altura, volume):
    #calcula a vazao maxima na barragem
    mmc = 0.0039 * volume ** 0.8122
    froe = 0.607 * (volume ** 0.295 * altura ** 1.24)

    return max(mmc, froe)

def qmax_secao(x, q_max_barr, volume):
    #calcula a vazao maxima em cada secao
    volume = volume * 10e-6
    if x == 0:
        return q_max_barr
    if volume > 6.2:
        return q_max_barr * 10 ** (-0.02/1609*x)
    else:
        a = 0.002 * np.log(volume) + 0.9626
        b = -0.20047 * (volume + 25000) ** -0.5979
        return q_max_barr * a * np.exp(b*x)

def pontos_tracado(linha, n=21):
    #definne n=21 pontos equidistantes na linha informada
    distances = np.linspace(0, linha.length, n)
    points = [linha.interpolate(distance) for distance in distances]

    return points, distances

def cotas(ponto_informado, srtm, altura):
    #calcula a cota de um ponto a partir do srtm
    cota_tal = list(srtm.sample([(ponto_informado.x, ponto_informado.y)]))
    cota_tal = cota_tal[0][0]
    cota_cor = cota_tal + altura

    return cota_tal, cota_cor

def simplificar_tracado(tracado, n):
    #simplifica o tracado do rio a partir de segmentos de retas
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
    #segmenta o tracado simplificado. necessario para o tracado das perpendiculares
    line_split = []
    coords = tracado.iloc[1]['geometry'].coords[:]
    for i in range(len(coords)):
        if i < len(coords)-1:
            a = coords[i]
            b = coords[i+1]
            l = LineString([a, b])
            line_split.append(l)

    return line_split

def perpendicular(linha, ponto, comprimento):
    #traca uma perpendicular a uma determinada reta em um determinado ponto
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
    #traca n=21 secoes perpendiculares de comprimento=400 metros equidistantes ao longo do tracado simplificado
    p, d = pontos_tracado(tracado.iloc[1]['geometry'], n)
    l = split_linha(tracado)
    tol = 1e-8
    for i, point in enumerate(p):
        for line in l:
            if line.distance(point) < tol:
                perp = perpendicular(line, point, comprimento)
                data = ['seçao {}'.format(i), d[i], perp]
                tracado.loc[len(tracado)] = data

    return tracado, d

def exportar_geopandas(tracado, nome_do_arquivo='tracado.shp'):
    #exporta o tracado como shp file
    tracado.crs = 'EPSG:31982'
    tracado = tracado.to_crs(epsg=4326)
    tracado.to_file(nome_do_arquivo)

def transformacao(x, y, d_to_m, new):
    #tranforma o sistema de coordenadas de um ponto de graus para metros ou vice versa
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
    #calcula a altimetria das secoes a partir do srtm
    cotas = []
    xs = []
    ys = []
    for linha in tracado.iloc[2:]['geometry']:
        p, d = pontos_tracado(linha, n=81)
        p = [i.coords[:][0] for i in p]
        x, y = [i[0] for i in p], [j[1] for j in p]
        xs.append(x)
        ys.append(y)
        pt = transformacao(x, y, d_to_m=False, new=True)
        pt = [(y, x) for x,y in zip(pt[0], pt[1])]
        cota = list(srtm.sample(pt))
        cota = [k[0] for k in cota]
        cotas.append(cota)

    return cotas, d, xs, ys

def line_coef(p1, p2):
    #calcula os coeficientes a, b de uma reta dado dois pontos
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    a = (y2 - y1) / (x2 - x1)
    b = y1 - a * x1     
    
    return a, b

def increase_resolution(x, y, n=100):
    #aumenta a resolução do perfil da secao
    x, y, = np.array(x), np.array(y)
    
    ixs = []
    iys = []

    for i in range(len(x)):
        if i < len(x)-1:
            ix = np.linspace(x[i], x[i+1], n)
            p1, p2 = (x[i], y[i]), (x[i+1], y[i+1])
            a, b = line_coef(p1, p2)
            iy = a*ix+b

            ixs = ixs + list(ix)
            iys = iys + list(iy)

    ixs, iys = np.array(ixs), np.array(iys)

    return ixs, iys

def raio_hidraulico(y, x, h_max):
    #calcula a area e o raio hidraulico de cada secao
    x, y = increase_resolution(x, y)
    
    yt = -1 * y + max(y)
    hs = np.linspace(0, h_max, 11) #11 alturas entre 0 e a altura maxima
    areas = []
    radius = []
    
    for h in hs[1:]: #para dez alturas entre um minumo e a altura maxima
        ytt = yt - (max(yt) - h)
        f = ytt > 0
        ytt, xt = ytt[f], x[f]
        
        #calcula a area pelo metodo dos trapezios
        area = np.trapz(y=ytt, x=xt)

        #calcula o perimetro molhado
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
    #equacao de manning
    a, r = np.array(a), np.array(r)
    q = k*a*r**(2/3)*j**(1/2)

    return q

def polyfit(x, y, x_i):
    #ajusta uma polinomial de terceiro grau
    coefs = poly.polyfit(x, y, 3)
    ffit = poly.polyval(x_i, coefs)

    return ffit

def clip_raster(secs, srtm, out_file):
    #corta o srtm a partir do tracado e das secoes
    secs.crs = 'EPSG:31982'
    secs = secs.to_crs(epsg=4326)
    minx, miny, maxx, maxy = min(secs.bounds['minx']), min(secs.bounds['miny']), max(secs.bounds['maxx']), max(secs.bounds['maxy'])
    bbox = box(minx, miny, maxx, maxy)
    out_image, out_transform = rasterio.mask.mask(srtm, [bbox], crop=True, nodata=-999)
    out_meta = srtm.meta
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})
    with rasterio.open(out_file, 'w', **out_meta) as dest:
        dest.write(out_image)

def get_coordinates(clipado):
    #retorna as coordenadas e a cota do srtm cortado
    w = clipado.width
    h = clipado.height
    ij = []
    for i in range(h):
        for j in range(w):
            a = clipado.xy(i, j, offset='center')
            ij.append(a)
    b = clipado.read(1).flatten()
    nanfilter = b > 0
    x, y = np.array([i[0] for i in ij]), np.array([j[1] for j in ij])
    ij = transformacao(y, x, d_to_m=True, new=True)
    xcoords = np.array(ij[0])
    ycoords = np.array(ij[1])

    return xcoords[nanfilter], ycoords[nanfilter], b[nanfilter]

def altura_de_agua_secoes(ds, dp, c, qmax_barr, v, h_barr):
    #calcula a altura de agua em cada secao
    ct = [i[40] for i in c]
    j = ((ct[0]+h_barr) - ct[-1])/ds[-1]

    #altura maxima simulada
    fc = 1 #fator de correcao (1 - 6) foi deifinido como 1 para todas as secoes para fins de simplificacao
    alt_max = h_barr/fc

    qs = []
    for i in ds:
        qs.append(qmax_secao(i, qmax_barr, v))

    areas = []
    raios = []
    for cotas in c:
        a, r, h = raio_hidraulico(cotas, dp, alt_max)
        areas.append(a)
        raios.append(r)
    alturas = h

    alturas_secoes = []
    for idx in range(len(areas)):
        qs_s = []
        for idx1 in range(len(areas[idx])):
            q = manning(areas[idx][idx1], raios[idx][idx1], j)
            qs_s.append(q)
        
        qs_s = np.insert(qs_s,0,0)
        a = polyfit(qs_s, alturas, qs[idx])
        a = a + ct[idx]
        alturas_secoes.append(a)

    return alturas_secoes, qs

def surfaces_to_kml(surf_surface, surf_water, flname):
    intersection, s1_split, s2_split = surf_surface.intersection(surf_water)
    linestrings = []
    number_of_previous_points=0

    for i in range(intersection.number_of_cells):
        index_to_find_length_of_line = i + number_of_previous_points
        number_of_points_of_line = intersection.lines[index_to_find_length_of_line]
        values = [intersection.lines[index_to_find_length_of_line+i+1] for i in range(number_of_points_of_line)]
        points = [intersection.points[value] for value in values]
        number_of_previous_points = number_of_previous_points + number_of_points_of_line
        linestrings.append(LineString(np.array(points)))

    multi_line = MultiLineString(linestrings)
    int_gdf = geopandas.GeoDataFrame(columns=['Nome', 'geometry'])
    int_gdf['Nome'] = ['Mancha de inundação']
    int_gdf['geometry'] = [multi_line]
    
    exportar_geopandas(int_gdf, flname)