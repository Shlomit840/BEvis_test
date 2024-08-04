#!/usr/bin/env python
# coding: utf-8

#--------------------------------------------

#impot configuration file

import yaml

config = yaml.safe_load(open("./BEvis_conf.yaml"))

# step_1: import packages and upload models.
#==============================================
# User input:

models_directory = config['step1']['models_directory']

#==============================================

import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

import networkx as nx

import os

# step_2: prepare data for COMETS

#==============================================

# User input:

nc = config['step2']['nc']

input_dir = config['step2']['input_dir']

cmd = "mkdir -p " + input_dir
os.system(cmd)

met_size_log_factor = config['step2']['met_size_log_factor']

met_size_over_all_factor = config['step2']['met_size_over_all_factor']

# step_4: Preparing  a dictionary for bacteria biomass, and compound amount, in the selected day, for using as nodes size.
#==============================================
# User input:

bm_factor = config['step4']['bm_factor']

#==============================================

def fix_cycle_near_10(n):

    #get day near 10, e.g., 79 to 80, 71 to 70

    #this is becuase days are in steps of 10 days

    t = n

    remider = n % 10

    if remider >= 5 :

        t = ((n // 10) +  1 )*10

    if remider >= 1  and remider < 5:
        t = ((n // 10))*10

    if n // 10 == 0 :

        t = 10

    return t

# for dictionary Biomass time

df_BM = pd.read_csv(input_dir + "df_BM.csv") 

#find day(cycle)

biomass_bool = df_BM.loc[:,"biomass"] == df_BM.loc[:,"biomass"].max()

#find max cycle (first index) with true biomass_bool

max_cycle  = df_BM.loc[biomass_bool, "cycle"].min()

# set day as nc cycles from max_cycle

day = fix_cycle_near_10(max_cycle - nc)

#day = 95

print("found day, ", day)

#added by hand to get results as example

df_BM_day = df_BM[df_BM.cycle ==day]

df_BM_day.set_index('species', inplace=True) 

def calculation_bm(row):

   cc = row['biomass']*int(float(bm_factor))

   if cc > 0:

      ccc = (np.log2(cc))*10

      if ccc > 0:

         return ccc

      else:

         return 1
   else:

      return 1

df_BM_day['bmcalc'] = df_BM_day.apply(calculation_bm, axis=1)

bac_bm_calc = df_BM_day.loc[:,'bmcalc'] 

bac_BM_dict = bac_bm_calc.to_dict() # dict of bacterium biomass representative value in time.

#=======================================

# for dictionary:   compound_id: [compound_name, formula] 

dicComp_df1 = pd.read_csv(input_dir + 'dicComp_df.csv')

dicComp_df2 = dicComp_df1.iloc[:,1:].T

dicComp_df2.columns = ["name", 'formula']

dicComp_df2['nameform'] = dicComp_df2.agg(list, axis=1) 

dic_namform = dict(dicComp_df2['nameform'])

#=================================================

# for dictionary metabolites amount in time

df_mets = pd.read_csv(input_dir + 'df_metabolites.csv')

if 'Unnamed: 0' in df_mets.columns:

    df_mets_ = df_mets.drop('Unnamed: 0',axis=1)

else:

    df_mets_ =  df_mets

name_cols = ["cycle"]

for col in df_mets_.columns[1:]:

    nam = dic_namform[col][0]

    name_cols.append(nam)

df_mets_.columns = name_cols

df_mets_day = df_mets_[df_mets.cycle ==day]

df_mets_day.set_index('cycle', inplace=True)

df_met_dayT = df_mets_day.T 

df_met_dayT.columns = ["vv"] #

def met_size(row):

   if row['vv'] > 0:

      cc = (np.log2(row['vv']*met_size_log_factor))*met_size_over_all_factor # להתאים <======

      if cc > 0:

         return cc

      else:

         return 1
   else:
      return 1
   
df_met_dayT['metcalc'] = df_met_dayT.apply(met_size, axis=1)

met_calc = df_met_dayT.loc[:,'metcalc'] 

met_calc_dict = met_calc.to_dict() # # dict of compounds representative value in time.

#=========================================================

#=========================================================

# for dictionary:   bacteria_name: [Exchanges_flux_table]

dic_EX_files = {} 

with open(input_dir + "model_id_list.txt") as file:

    lines = [line.rstrip() for line in file]

for file_name in lines:

    mod = pd.read_csv(input_dir + "/model_exchanges/"  + file_name)

    mod = mod.drop(columns = ['Unnamed: 0'])

    mod.set_index('cycle', inplace=True)

    mod1 = mod.loc[:, (mod != 0).any(axis=0)] 

    C_cols = [] # just compounds with carbone.

    for col in mod1.columns:

        col2 = col.replace('EX_','')

        if 'cpd11416' not in col:

            if col2 in dic_namform.keys():

                formula = dic_namform[col2][1]

                if formula[0] == 'C' and formula[1].islower() == False:

                    C_cols.append(col)   

    mod3 = mod1[C_cols]

    bac_name = file_name

    bac1_name__ = bac_name.replace('.csv','')

    dic_EX_files[bac1_name__] = mod3 

#===================================

#===================================

# for one cycle, tuple lists of (bacteria_name,exchange_compound_name,flux).

# devide between compounds uptaked and compouns secrated.

#day = 90

secretions = []

uptakes = []

for k in dic_EX_files.keys():

    df_flux = dic_EX_files[k]

    loc = df_flux.index.get_loc(day)

    df_flux_day = df_flux.iloc[loc,:]

    df_flux_day_ = pd.DataFrame(df_flux_day)

    df_flux_day_.columns = ['flux']

    df_EX_secretion = df_flux_day_[(df_flux_day_["flux"] > 0.0015)]

    df_EX_uptake = df_flux_day_[(df_flux_day_['flux'] < -0.0015)]

    met_sec = df_EX_secretion.index

    for i in met_sec:

        tup = (k,i,float(df_EX_secretion.loc[i]))

        secretions.append(tup)

    met_uptake = df_EX_uptake.index

    for i in met_uptake:

        tup = (i, k, float(df_EX_uptake.loc[i]))

        uptakes.append(tup)

#==============================

secretionsN = [] 

for tup in secretions:

    id = tup[1].replace('EX_','')

    name = dic_namform[id][0]

    tupN = (tup[0], name, tup[2])

    secretionsN.append(tupN)      
    
uptakesN = [] 

for tup in uptakes:

    id = tup[0].replace('EX_','')

    name = dic_namform[id][0]

    tupN = (name, tup[1], tup[2])

    uptakesN.append(tupN)
         
#==================================

#day #from end of step 2

output_dir = config['step3']['output_dir']
output_dir
cmd = "mkdir -p " + output_dir

os.system(cmd)

uptake_of_interest = config['step3']['uptake_of_interest']

secretions_of_interest = config['step3']['secretions_of_interest']

# Pick tuples for the figure.

uptake_of_interest_tup = []

for tup in uptakesN:

    for b in uptake_of_interest:

        if b in tup:

            uptake_of_interest_tup.append(tup)

secretions_of_interest_tup = []

for tup in secretionsN:

    for i in secretions_of_interest:

        if i in tup:

            secretions_of_interest_tup.append(tup)

secretion_in_uptake = []

for tups in secretionsN:

    for tupu in uptakesN:

        tupsl = list(tups)

        if tupsl[1] in tupu:

            if tups not in secretion_in_uptake:

                secretion_in_uptake.append(tups)

uptake_in_secretion = []

for tups in secretionsN:

    for tupu in uptakesN:

        tupul = list(tupu)

        if tupul[0] in tups:

            if tupu not in uptake_in_secretion:

                uptake_in_secretion.append(tupu)

#===========================================

# Preparing a dictionary of bacterium and compound as key, and flux as value.

dic = {} 

for tup in secretion_in_uptake:

    ltup = list(tup)

    lltup = ltup[0:2]

    key_secretion = tuple(lltup)

    dic[key_secretion] = ltup[2]

for tup2 in uptake_in_secretion:
    
    ltup2 = list(tup2)

    lltup2 = ltup2[0:2]

    key_secretion = tuple(lltup2)

    dic[key_secretion] = ltup2[2]

for tup3 in secretions_of_interest_tup:

    ltup = list(tup3)

    lltup = ltup[0:2]

    key_secretion = tuple(lltup)

    dic[key_secretion] = ltup[2]

for tup4 in uptake_of_interest_tup:

    ltup2 = list(tup4)

    lltup2 = ltup2[0:2]

    key_secretion = tuple(lltup2)

    dic[key_secretion] = ltup2[2] 


#================================================== 

# step_5: define edges width, define vertical groups (layers) of bacteria and metabolites, and define color and size of nodes.

#==============================================

# User input:

edges_factor_to_devide = config['step5']['edges_factor_to_devide']

# vertical groups (layers)

# BTEX_degradors # layer 1

bac_A = config['step5']['bac_A']

# Fermentors layer 3

bac_B = config['step5']['bac_B']

# Methanogens # layer 5

bac_C = config['step5']['bac_C']

#==============================================


G = nx.DiGraph()

G.add_edges_from(dic.keys())

# defining edges width by flux.

wid = []

for edge in G.edges:

    wid.append(dic[edge])

big_thick = []

for x in wid:

    xx = abs(x) * 1000

    big_thick.append(xx)

log_thick = []

for n in big_thick:

    if n > 0:

        log_thick.append(np.log2(n))

    else:

        log_thick.append(0.1)

log_thick_2 = []

for n in log_thick:

    log_thick_2.append(n/edges_factor_to_devide)

# compounds for vertical groups (layers)

comp_B = [] # layer 2

for item in dic.items():

    if item[0][0] in bac_A:

        comp_B.append(item[0][1])

comp_C = [] # layer 2

for item in dic.items():

    if item[0][1] in bac_B:

        comp_C.append(item[0][0])

comp_D = [] # layer 4

for item in dic.items():

    if item[0][0] in bac_B:

        comp_D.append(item[0][1])

comp_E = [] # layer 4

for item in dic.items():

    if item[0][1] in bac_C:

        comp_E.append(item[0][0])

comp_F = [] # layer 6

for item in dic.items():

    if item[0][0] in bac_C:

        comp_F.append(item[0][1])

#=====================

# vertical groups (layers)

nodes_from_layers = []

for i in uptake_of_interest:

    if i in list(G.nodes):

        G.add_node(i,layer=0)

        nodes_from_layers.append(i)

for i in bac_A: 

    if i in list(G.nodes):

        G.add_node(i,layer=1)

        nodes_from_layers.append(i)

for i in comp_B:

    if i in list(G.nodes):

        if i not in nodes_from_layers:

            G.add_node(i,layer=2)

            nodes_from_layers.append(i)

for i in comp_C:

    if i in list(G.nodes):

        if i not in nodes_from_layers:

            G.add_node(i,layer=2)

            nodes_from_layers.append(i)

for i in bac_B: 

    if i in list(G.nodes):

        G.add_node(i,layer=3)

        nodes_from_layers.append(i)

for i in comp_D:

    if i in list(G.nodes):

        if i not in nodes_from_layers:

            G.add_node(i,layer=4)

            nodes_from_layers.append(i)

for i in comp_E:

    if i in list(G.nodes):

        if i not in nodes_from_layers:

            G.add_node(i,layer=4)

            nodes_from_layers.append(i)

for i in bac_C: 

    if i in list(G.nodes):

        G.add_node(i,layer=5)

        nodes_from_layers.append(i)

for i in comp_F:

    if i in list(G.nodes):

        if i not in nodes_from_layers:

            G.add_node(i,layer=6)

            nodes_from_layers.append(i)

for i in list(G.nodes):

    if i not in nodes_from_layers:

        G.add_node(i,layer=7)

#==================================

# defining nodes color and size.

color_map = []

node_sizes = []

for node in G.nodes:

    if node in bac_A:

        color_map.append('mediumaquamarine')
        
        node_sizes.append(bac_BM_dict[node])

    elif node in  bac_B:

        color_map.append('burlywood')

        node_sizes.append(bac_BM_dict[node])

    elif node in bac_C:

        color_map.append('darkorange')

        node_sizes.append(bac_BM_dict[node])

    elif node in uptake_of_interest:

        color_map.append('brown')

        node_sizes.append(met_calc_dict[node])

    elif node in secretions_of_interest:

        color_map.append('red')

        try:

            node_sizes.append(met_calc_dict[node])

        except:

            node_sizes.append(10)
    else:

        try:

            node_sizes.append(met_calc_dict[node])

            color_map.append('deepskyblue')

        except:

            node_sizes.append(10)

            color_map.append('grey')
           
#=========================================

# step_6: create a plot

#==============================================

# User input:

plot_filename = config['step6']['plot_filename']

#==============================================

pos = nx.multipartite_layout(G, subset_key="layer")

plt.figure(figsize=(20, 15))

edges = nx.draw_networkx_edges(G, pos,width= log_thick_2, alpha=0.5)

nodes = nx.draw_networkx_nodes(G, pos, node_color=color_map, node_size=node_sizes, alpha=1)

labels = nx.draw_networkx_labels(G, pos, font_size=10)

plt.savefig(output_dir + plot_filename, dpi=600) # Change file name.
