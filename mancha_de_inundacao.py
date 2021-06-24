#!/usr/bin/env python

#importando os pacotes
import geopandas
import matplotlib.pyplot as plt

#definindo caminho dos arquivos necessários
bacias_path = 'deps/bacias/Bacia_Hidrografica.shp'
areas_urbanizadas_min = 'deps/areas_min_rs/areas_min_rs.shp'
areas_urbanizadas_max = 'deps/areas_max_rs/areas_max_rs.shp'

fig, ax = plt.subplots(figsize=(8,8))

gp_bacias = geopandas.read_file(bacias_path)
gp_bacias.plot(ax=ax, color='gray')

gp_amin = geopandas.read_file(areas_urbanizadas_min)
gp_amin.plot(ax=ax, color='green')

gp_amax = geopandas.read_file(areas_urbanizadas_max)
gp_amax.plot(ax=ax, color='red')

plt.show()