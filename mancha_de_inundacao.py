"""mancha_de_inundacao.py: modulo com as funcoes necessarias para o programa mancha_de_inundacao."""

__author__ = "Roberto Mentzingen Rolo"

#importando pacotes
import numpy as np
import numpy.polynomial.polynomial as poly
from shapely.geometry import Point, LineString, box
from shapely import affinity
from pyproj import Transformer, transform
import math
import rasterio.mask
import rasterio
from scipy.interpolate import Rbf
import simplekml

def rotate_l(l1, drange):
    #gira um linha em seu centro um valor aleatorio em graus entre um range definido
    g = np.random.uniform(drange[0],drange[1])
    nl = affinity.rotate(l1, g, 'center')
    return nl

def rotate_secs(sec_df, maxiter=5000, drange=[-5,5]):
    #rotaciona as secoes ate que nao haja mais intersecoes 
    ndf = secs.iloc[2:].copy(deep=False)
    ndf = ndf.reset_index(drop=True)
    inter = True
    
    while inter and maxiter > 0:
        maxiter = maxiter -1
        inter = False

        for idx in range(ndf.shape[0]-1):
            if idx == 0:
                l1 = ndf.iloc[idx].geometry
                l2 = ndf.iloc[idx+1].geometry
                if l1.intersects(l2):
                    nl = rotate_l(l1, drange)
                    ndf.iloc[idx].geometry = nl
                    inter=True
            elif idx == ndf.shape[0]-1:
                l1 = ndf.iloc[idx].geometry
                l2 = ndf.iloc[idx-1].geometry
                if l1.intersects(l2):
                    nl = rotate_l(l1, drange)
                    ndf.geometry.iloc[idx] = nl
                    inter=True
            else:
                for i in [-1,1]:
                        l1 = ndf.iloc[idx].geometry
                        l2 = ndf.iloc[idx+i].geometry
                        if l1.intersects(l2):
                            nl = rotate_l(l1, drange)
                            ndf.geometry.iloc[idx] = nl
                            inter=True
                    
    print('Iteração: ', maxiter)
    print('interseção: ', inter)
    return ndf

def crio(volume):
    #calcula o comprimento do rio a ser modelado a partir da barragem
    volume = volume * 10e-6 # tranformacao para hm
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

def rbf_interpolation(x, y, v, xi, yi, function='linear'):
    #interpola a cota das secoes a partir de uma funcao linear
    x, y, z, d = x, y, np.zeros(len(x)), v
    rbfi = Rbf(x, y, z, d, function=function)
    di = rbfi(xi, yi, np.zeros(len(xi)))

    return di

def points_to_kml(x, y, mancha, flname):
    #salva o conjunto de pontos em formato kml para visualizacao no google earth
    x = np.array(x)
    y = np.array(y)
    f = mancha == 1
    x = x[f]
    y = y[f]
    xy = transformacao(x, y, d_to_m=False, new=True)

    kml = simplekml.Kml()
    for x, y in zip(xy[0], xy[1]):
        pnt = kml.newpoint(description='ponto inundado', coords=[(y, x)])
        pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/water.png'
        pnt.style.iconstyle.scale = 0.5

    kml.save(flname)
