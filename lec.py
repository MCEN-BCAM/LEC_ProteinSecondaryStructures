'''
	Author: Rodrigo A. Moreira (C) 2023
	https://orcid.org/0000-0002-7605-8722
	LICENSE: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/)
'''

import math
import networkx as nx

def read_network(net):
	G = nx.Graph()
	with open(net,'r') as file:
		line = file.readline()
		G.add_nodes_from(range(1,int(line)+1))
		line = file.readline()
		while line:
			splt = line.split()
			d = math.floor(float(splt[2])*10 + 0.5)/10.0
			G.add_edge( int(splt[0]),int(splt[1]), weight=d )
			line = file.readline()
	return G

def kappa(net,max):
	'''
	INPUT: numbers of nodes and a list of edges
	OUTPUT: Partition of Euler Characteristic for each node in the clique complex
	'''
	def filtration():
		G = net.copy()
		G.remove_edges_from([(n1, n2) for n1,n2,w in G.edges(data="weight") if w > max])
		return G
	KAPPA = [1.0 for k in net.nodes()]
	for clq in nx.enumerate_all_cliques(filtration()):
		lclq = len(clq)
		if lclq > 1:
			ilclq = 1.0/lclq
			sgn = -1 if lclq%2 == 0 else 1
			for k in clq:
				KAPPA[k-1] += sgn*ilclq
	return KAPPA

def listtostr(l):
	str = ''
	for k in l:
		str += ' %.3f'%k
	return str+"\n"

import sys
PDB = sys.argv[1]

print("Processing ",PDB)
NET = read_network('%s.pdb.network_backboneRE_heavy_gt2'%PDB)
L = ['0.000']
for n1,n2,w in NET.edges(data='weight'):
	k = '%.1f'%w
	if k not in L:
		L.append(k)
L = [float(k) for k in L]
list.sort(L)
pfile1 = open('%s.kappas'%PDB,'w')
pfile2 = open('%s.reslec'%PDB,'w')
for max in L:
	k = kappa(NET,max)
	pfile1.write('%.3f'%max+listtostr(k))
	reslec = [k[i]+k[i+1]+k[i+2] for i in range(0,len(k),3)]
	pfile2.write('%.3f'%max+listtostr(reslec))
pfile1.close()
pfile2.close()

