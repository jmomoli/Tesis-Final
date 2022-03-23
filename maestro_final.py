### Librerias ###
import readdy
import sys
import numpy as np
import os
import csv
import random

### Definiciones ###
print("Instrucciones: Para ejecutar este codigo en la consola, debe agregar al comando el numero de celulosomas que desea simular y el tipo del mismo(BLS/NATURAL/LIBRES), y que enzimas (endo, exo, endo_exo)")

n = int(sys.argv[1]) #Numero de celulosomas o enzimas libres, insertar un numero mayor a 0
clase = sys.argv[2] #Tipo de celulosomas(BLS/NATURAL/LIBRES)
enzimeus = sys.argv[3]
### Caja ###
system = readdy.ReactionDiffusionSystem(box_size=[60, 60, 60], periodic_boundary_conditions=[True, True, True])

### Funciones ###
#Generales#
def diccionario_particulas(l,r):
	"""Recible una lista de particulas y crea el diccionario con dichas particulas y sus radios para poder ser visualizados luego en el VMD"""
	diccionario = {}
	for i in range(len(l)):
		diccionario[l[i]] = r
	return diccionario

def add_topology_species(l,d):
	"""La funcion recibe una lista de strings correspondiente a la cantidad de topologias a agregar (l) y le asigna un valor de constante de difussion (d)"""
	#system.add_topology_species("BLS1", diffusion_constant = d)
	for i in range(len(l)):# el if sirve para cuando tenemos los primeros monomeros reemplzados por BLS1 para calcular el MSD
            system.add_topology_species(l[i], diffusion_constant = d)
        #	if l[i] == "BLS1":
	#		print("LORD")
	#	else:
	#		system.add_topology_species(l[i], diffusion_constant = d)

#BLS#
def listar_BLS(n):
	"""La funcion recibe un numero n de celulosomas y crea una lista con 10 strings desde BLS1 a BLS(n*10)."""
	BLS = []
	for i in range(1,(n*10+1)):
		BLS.append("BLS" + str(i))
	return BLS

def harmonic_bond_BLS(l,f=10):
	"""Recibe una lista l de strings correspondiente a la cantidad de particulas de BLS y crea los enlaces armonicos correspondientes con fuerza f"""
	for i in range(0,len(l),10):
	#Saltea de 10 en 10 para ir de celulosoma en celulosoma#
		BLS = l[i:i+10] #Creo un slice de las 10 subunidades del celulosoma, esto me permite trabajar con un unico celulosoma a la vez#
		system.topologies.configure_harmonic_bond(BLS[0], BLS[4], force_constant = f, length = (BLS1 + BLS1)) #Cierro el primero anillo#
		system.topologies.configure_harmonic_bond(BLS[5], BLS[9], force_constant = f, length = (BLS1 + BLS1)) #Cierro el segundo anillo#
		for j in range(len(BLS)-1):
		#Itera en la lista BLS#
			system.topologies.configure_harmonic_bond(BLS[j], BLS[j+1], force_constant = f, length = (BLS1 + BLS1)) #Crea en enlace entre cada unidad de BLS y su siguiente#
		for k in range(0,5):
		#Itero en en primer anillo#
			system.topologies.configure_harmonic_bond(BLS[k], BLS[k+5], force_constant = f, length = (BLS1 + BLS1)) #Crea en enlace entre cada subunidad del primer anillo con su contraparte del segundo anillo#

def harmonic_angle_BLS(l,f=100,e1=1.8849,e2=1.57079):
	"""Similar a harmonic_bonds pero esta crea los angulos simples entre los monomeros de BLS.
	La fuerza esta determinada por f y los angulos de equilibrios para un mismo anillo y para con el otro anillo estan determinados por e1 y e2 respectivamente.
	Los subindices pueden parecer un lio pero esta checkeado que se crean los angulos correspondientes."""
	for i in range(0,len(l),10):
		BLS = l[i:i+10]
		system.topologies.configure_harmonic_angle(BLS[3],BLS[4],BLS[0], force_constant = f, equilibrium_angle = e1)
		system.topologies.configure_harmonic_angle(BLS[4],BLS[0],BLS[1], force_constant = f, equilibrium_angle = e1)
		system.topologies.configure_harmonic_angle(BLS[9],BLS[5],BLS[6], force_constant = f, equilibrium_angle = e1)
		system.topologies.configure_harmonic_angle(BLS[8],BLS[9],BLS[0], force_constant= f, equilibrium_angle = e1)
		system.topologies.configure_harmonic_angle(BLS[9],BLS[5],BLS[6], force_constant= f, equilibrium_angle = e1)
		system.topologies.configure_harmonic_angle(BLS[8],BLS[9],BLS[4], force_constant = f, equilibrium_angle = e2)
		system.topologies.configure_harmonic_angle(BLS[3],BLS[4],BLS[9], force_constant = f, equilibrium_angle = e2)
		system.topologies.configure_harmonic_angle(BLS[4],BLS[0],BLS[5], force_constant = f, equilibrium_angle = e2)
		for j in range(0,3):
			system.topologies.configure_harmonic_angle(BLS[j], BLS[j+1], BLS[j+2], force_constant = f, equilibrium_angle = e1)
			system.topologies.configure_harmonic_angle(BLS[j+5], BLS[j+6], BLS[j+7], force_constant = f, equilibrium_angle = e1)
			system.topologies.configure_harmonic_angle(BLS[j], BLS[j+1], BLS[j+6], force_constant = f, equilibrium_angle = e2)
			system.topologies.configure_harmonic_angle(BLS[j+5], BLS[j+6], BLS[j+1], force_constant = f, equilibrium_angle = e2)
def dihedral_BLS(l,f=100):
	""" Setea los angulos dihedricos a partir de una lista l y una fuerza f. Los subindices pueden parecer un lio pero esta checkeado que se crean los angulos correspondientes."""
	for i in range(0,len(l),10):
		BLS = l[i:i+10]
		system.topologies.configure_cosine_dihedral(BLS[2],BLS[3],BLS[4],BLS[0], force_constant=f,multiplicity=1.,phi0=0.)
		system.topologies.configure_cosine_dihedral(BLS[3],BLS[4],BLS[0],BLS[1], force_constant=f,multiplicity=1.,phi0=0.)
		system.topologies.configure_cosine_dihedral(BLS[4],BLS[0],BLS[1],BLS[2], force_constant=f,multiplicity=1.,phi0=0.)
		system.topologies.configure_cosine_dihedral(BLS[7],BLS[8],BLS[9],BLS[5], force_constant=f,multiplicity=1.,phi0=0.)
		system.topologies.configure_cosine_dihedral(BLS[8],BLS[9],BLS[5],BLS[6], force_constant=f,multiplicity=1.,phi0=0.)
		system.topologies.configure_cosine_dihedral(BLS[9],BLS[5],BLS[6],BLS[7], force_constant=f,multiplicity=1.,phi0=0.)
		for j in range(0,2):
			system.topologies.configure_cosine_dihedral(BLS[j],BLS[j+1],BLS[j+2],BLS[j+3], force_constant=f,multiplicity=1.,phi0=0.)
			system.topologies.configure_cosine_dihedral(BLS[j+5],BLS[j+6],BLS[j+7],BLS[j+8], force_constant=f,multiplicity=1.,phi0=0.)
def repulsivos_BLS(l,f = 10):
	"""Recorre una lista seteando el potencial repulsivo entre todas las particulas de la lista.
	Es el volumen de exclusion."""
	for i in range(len(l)+1):
		for j in range(i+1,len(l)):
			system.potentials.add_harmonic_repulsion(l[i],l[j], force_constant = f, interaction_distance = BLS1+BLS1)

#Natural#
def listar_scaffold(n):
	"""La funcion recibe un numero n de celulosomas y crea una lista con 10 strings desde scaffold1 a scaffold(n*10)."""
	scaffold = []
	for i in range(1,(n*10+1)):
		scaffold.append("scaffold" + str(i))
	return scaffold
def harmonic_bond_scaffold(l,f=5):
	"""Similar a su analogo para BLS, genera los enlaces del celulosoma natural recibiendo una lista de strings y un valor de fuerza de union."""
	for i in range(0,len(l),10):
		scaffold = l[i:i+10]
		for j in range(0,len(scaffold)-1):
			system.topologies.configure_harmonic_bond(scaffold[j],scaffold[j+1], force_constant = f, length = (scaffold1 + scaffold1))
def harmonic_angle_scaffold(l, f=10, e = 3.14):
	"""Similar a su analogo para BLS, genera los algulos del celulosoma natural recibiendo una lista de strings, un valor de fuerza y un angulo"""
	for i in range(0,len(l),10):
		scaffold = l[i:i+10]
		for j in range(0,len(scaffold)-2):
			system.topologies.configure_harmonic_angle(scaffold[j],scaffold[j+1],scaffold[j+2], force_constant = f, equilibrium_angle = e)

#Linkers#
def listar_linkers(n):
	"""La funcion recibe un numero n de celulosomas y crea una lista con 10 strings desde F1 a F(n*10)."""
	F = []
	for i in range(1,(n*10+1)):
		F.append("F" + str(i))
	return F

def angle_linker(l1,l2,f=100,e = 1.8849):
	"""Crea los angulos entre los monomeros (BLS o naturales) y los linkers"""
	for i in range(0,len(l1),10):
		BLS = l1[i:i+10]
		F = l2[i:i+10]
		for j in range(5):
			system.topologies.configure_harmonic_angle(BLS[j],BLS[j+5],F[j+5], force_constant=f, equilibrium_angle = e)
			system.topologies.configure_harmonic_angle(BLS[j+5],BLS[j],F[j], force_constant=f, equilibrium_angle = e)

def BLS_linker(l1,l2,f=10):
	"""Crea los enlaces entre monomeros de BLS y los linkers. El orden de las listas es irrelevante"""
	for i in range(0,len(l1)):
		system.potentials.add_weak_interaction_piecewise_harmonic(l1[i],l2[i], force_constant= 100,desired_distance=1.3, depth=300, cutoff=5)
		#system.potentials.add_lennard_jones(l1[i],l2[i], m=4, n=2, cutoff=10, shift=True, epsilon=40, sigma=1.3)
		for j in range(0,len(l2)):
			if j != i:
				#system.potentials.add_harmonic_repulsion(l1[i], l2[j], force_constant=10, interaction_distance=(F1 + BLS1))
			#else:
				system.potentials.add_harmonic_repulsion(l1[i], l2[j], force_constant=50, interaction_distance=(F1 + BLS1))


def linker_enzima(l):
	for i in range(0,len(l)):
		for j in range(0,len(list_enzimas)):
				system.topologies.configure_harmonic_bond(l[i],list_enzimas[j], force_constant = 10, length =(F1+endo))
				system.potentials.add_harmonic_repulsion(l[i],list_enzimas[j],force_constant = 20, interaction_distance = (F1+endo))

def BLS_enzima(l):
	for i in range(0,len(l)):
		for j in range(0,len(list_enzimas)):
			system.potentials.add_harmonic_repulsion(l[i],list_enzimas[j],force_constant=50, interaction_distance = (endo+BLS1+0.3))



def BLS():
	"""Ejecuta las funciones necesarias para crear n decameros de BLS, con la salvedad la BLS_list que debe ser ejecutada con anterioridad para poder definir los radios"""
	print(BLS_list)	
	add_topology_species(BLS_list,DIFUSION)#Modificar aca la difusion
	harmonic_bond_BLS(BLS_list,5)
	harmonic_angle_BLS(BLS_list,100)
	dihedral_BLS(BLS_list,100)
	repulsivos_BLS(BLS_list,2.5)
	BLS_linker(BLS_list, F_list)
	angle_linker(BLS_list,F_list)
	linker_enzima(F_list)
	BLS_enzima(BLS_list)

def LIBRES():
	linker_enzima(F_list)

def main_config(clase = sys.argv[2]):
	tipo = clase
	if tipo == "BLS":
		print("Usted esta por simular "+sys.argv[1]+" celulosomas del tipo BLS")
		BLS()
	elif tipo == "NATURAL":
		print("Falta hacer funcion")
	elif tipo == "LIBRES":
		LIBRES()
	else:
		print("Error en el imput, los tipos validos son BLS/NATURAL/LIBRES")

def reemplazar(l):
	"""Esta funcion esta para reemplazar el primer monomero de cada decamero por BLS1 para calcular mas facil el MSD"""
	system.add_topology_species("BLS1", diffusion_constant = DIFUSION)
	for i in range(0,len(l),10):
		l[i] = "BLS1"

def difusion_set(clase):
    if clase =="BLS":
        return 5
    else:
        return 0.24

### Ejecucion ####
"""Crea las listas y asigna los radios necesarios para poder ejecutar main. Los linkers se listan siempre"""
DIFUSION = difusion_set(sys.argv[2])

if n<=0:
	print("Inserte un numero mayor a 0")
else:
	diccionario =  {"headA": 1, "A": 1, "core": 1, "coreB": 1, "B": 1, "headB": 1, "endo": 1, "endoB": 1, "cbm": 1, "cbmB": 1, "exo": 1, "exoB": 1, "exoI": 1, "bg": 0.5, "bgB": 0.5, "J": 0.5, "celobiosa": 1, "celobiosaB": 1,"celobiosaI": 1, "celobiosaR": 1, "celobiosaX": 1, "glucosa": 1}
	if clase=="BLS":
		BLS_list = listar_BLS(n)
		#reemplazar(BLS_list)
		diccionario.update(diccionario_particulas(BLS_list, 1))
		for i in range(1,len(BLS_list)+1):
			exec("BLS"+str(i) + "=" + "1.")
	if clase=="NATURAL":
		scaffold_list = listar_scaffold(n)
		diccionario.update(diccionario_particulas(scaffold_list,1))
		for i in range(1,len(scaffold_list)+1):
			exec("scaffold"+str(i) + "=" + "1")

	F_list = listar_linkers(n)
	add_topology_species(F_list,DIFUSION)
	diccionario.update(diccionario_particulas(F_list,0.3))
	for i in range(1,len(F_list)+1):
		exec("F"+str(i) + "=" + "0.3")

######## Topologies species ########

## cellulose ##
## core ##
system.add_topology_species("core", 0.)
system.add_topology_species("coreB", 0.)
## extremes ##
system.add_topology_species("headA", 0.)
system.add_topology_species("headB", 0.)
system.add_topology_species("A", 0.)
system.add_topology_species("B", 0.)
## cellobiose ##
system.add_topology_species("celobiosa", 1.)  # cellobiose
system.add_topology_species("celobiosaB", 1.)  # cellobiose bound to bg
system.add_topology_species("celobiosaI", 1.)  # cellobiose that inhibit exo

## enzymes ###
## endoglucanase ##
system.add_topology_species("endo", diffusion_constant=DIFUSION)  # enzyme not bound to cellulose
system.add_topology_species("endoB", diffusion_constant=DIFUSION)  # enzyme bound to cellulose
## exoglucanase ##
system.add_topology_species("exo", diffusion_constant=DIFUSION)  # enzyme not bound to cellulose
system.add_topology_species("exoB", diffusion_constant=DIFUSION)  # enzyme bound to cellulose
system.add_topology_species("exoI", diffusion_constant=DIFUSION)  # enzyme inhibited by product
## beta glucosidase ##
system.add_topology_species("bg", diffusion_constant=1.)  # enzyme not bound to cellobiose
system.add_topology_species("bgB", diffusion_constant=1.)  # enzyme not bound to cellobiose
## dummy enzyme ##
system.add_topology_species("J", diffusion_constant=1.)  # enzyme not bound to cellulose
## cellulose binding module ##
system.add_topology_species("cbm", diffusion_constant=1.)  # enzyme not bound to cellulose
system.add_topology_species("cbmB", diffusion_constant=1.)  # enzyme bound to cellulose

######## Not topologies species ########

## products ##
system.add_species("celobiosaX", diffusion_constant=1.)
system.add_species("celobiosaR", diffusion_constant=1.)
system.add_species("glucosa", diffusion_constant=1.)

################################ Adding Topologies ##############################

system.topologies.add_type("TopologyA")  # For scaffolds
system.topologies.add_type("TopologyB")  # For cellulose with no bound enzymes
system.topologies.add_type("TopologyC")  # For enzymes and cellulose with bound enzymes

################################ Particles radii ##############################

## cellulose ##
headA = 1.
headB = 1.
core = 1.
coreB = 1.
A = 1.
B = 1.

## enzymes ##
endo = 1.
endoB = 1.
exo = 1.
exoB = 1.
exoI = 1.
bg = 1.
bgB = 1.
cbm = 1.
cbmB = 1.
J = 1.

list_enzimas = ['endo','endoB','exo','exoB','exoI','bg','bgB','cbm','cbmB','J']

## products ##
celobiosa = 1.0
celobiosaI = 1.0
celobiosaB = 1.0
celobiosaR = 1.0
celobiosaX = 1.0
glucosa = 1.0

################################ Adding Potentials ##############################


## Atractive potentials ##


## Combinatoria de potenciales entre todas las especies descriptas a continuacion

## En el diccionario particulas los keys son los nombres de las especies y los values son los valores de los radios definidios

particulas = {"headA": 1, "A": 1, "core": 1, "coreB": 1, "B": 1, "headB": 1, "endo": 1, "endoB": 1, "cbm": 1, "cbmB": 1,
			  "exo": 1, "exoB": 1, "exoI": 1, "bg": 1, "bgB": 1, "J": 1, "celobiosa": 1, "celobiosaB": 1,
			  "celobiosaI": 1}

for k1 in range(0, len(particulas)):
	for k2 in range(k1, len(particulas)):
		system.topologies.configure_harmonic_bond(list(particulas)[k1], list(particulas)[k2], 10.,
												  (particulas[list(particulas)[k1]] + particulas[list(particulas)[k2]]))

## Repulsive potentials ##


## Combinatoria de potenciales entre todas las especies descriptas a continuacion

##### MAIN MAIN MAIN ###
main_config()



for k1 in range(0, len(diccionario)):
            for k2 in range(k1, len(diccionario)):
                    system.potentials.add_harmonic_repulsion(list(diccionario)[k1], list(diccionario)[k2], 10., (
                                            diccionario[list(diccionario)[k1]] + diccionario[
                                    list(diccionario)[k2]]))  # esta en 1000 ERA 																									500 respecto sistema





################################ Adding reactions ##############################

## exoglucanase ##

system.topologies.add_spatial_reaction('exo bind A: TopologyB(headA) + TopologyC(exo) -> TopologyC(A--exoB)', rate=0.1,
									   radius=exo + headA)  # (k1) binding of exo to headA (cellulose topology B)

system.topologies.add_spatial_reaction('exo bind B: TopologyB(headB) + TopologyC(exo) -> TopologyC(B--exoB)', rate=0.1,
									   radius=exo + headA)  # (k1) binding of exo to headB (cellulose topology B)

system.topologies.add_spatial_reaction('exo bind AA: TopologyC(headA) + TopologyC(exo) -> TopologyC(A--exoB)', rate=0.1,
									   radius=exo + headA)  # (k1) binding of exo to headA (cellulose topology C)

system.topologies.add_spatial_reaction('exo bind BB: TopologyC(headB) + TopologyC(exo) -> TopologyC(B--exoB)', rate=0.1,
									   radius=exo + headA)  # (k1) binding of exo to headB (cellulose topology C)


def exog1_rate_function(topology):  # kcat para reaccion de exoglucanasa
	return 0.01


def exog1_reaction_function(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertices = topology.get_graph().get_vertices()
	edges = topology.get_graph().get_edges()
	lista_core = []  # lista con index de core en orden
	lista_coreB_index = []  # lista con index de coreB en orden
	lista_coreB = []
	lista_endoB = []  # lista con index de endoB en orden
	lista_endoB_index = []
	lista_endo = []
	lista_exoB = []  # lista con index de endoB en orden
	lista_exo = []
	lista_core_coreB = []  # lista con index de core y coreB en orden
	lista_coreB_endoB = {}  # diccionario con keys = nombre de particula, values = listas con index de particulas en orden
	lista_headB = []  # lista con index de headB en orden
	lista_headA = []
	lista_A_B = []
	lista_A_B_index = []
	lista_A = []
	lista_B = []
	lista_A_B_exoB = []
	count = 0
	for v in vertices:
		a = topology.particle_type_of_vertex(v)
		print(a, count, topology.particle_id_of_vertex(v))
		if a == 'headB':
			lista_headB.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headB es", lista_headB)
			largocadena = count + 1  # calculo de largo de cadena. Sirve para limitar la reaccion a cierto largo de cadena
			anteultimo = largocadena - 2  # calculo de anteultimo valor de cadena. Sirve para que endo no reaccione con borde generando glucosa.
			# print("largo cadena es " + str(largocadena))
			count = count + 1
		elif a == 'headA':
			lista_headA.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headA es", lista_headA)
			count = count + 1
		elif a == 'A':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_A_B.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			lista_A.append(count)
			lista_A_B_index.append(count)
			# lista_A_index.append(count)
			# lista_A_B_exoB[a0] = lista_A_B		#agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'B':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_A_B.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			lista_B.append(count)
			lista_A_B_index.append(count)
			# lista_A_index.append(count)
			# lista_core_coreB.append(count)
			# lista_A_B_exoB[a0] = lista_A_B		#agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'core':
			lista_core.append(count)  # agrega valores de index de las particulas core a lista_core
			lista_core_coreB.append(count)
			# print("lista core es ", lista_core)
			count = count + 1
		elif a == 'coreB':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_coreB.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			lista_coreB_index.append(count)
			lista_core_coreB.append(count)
			lista_coreB_endoB[
				a0] = lista_coreB  # agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'exo':
			lista_exo.append(count)
			count = count + 1
		elif a == 'exoB':
			# print(v)
			a1 = str(a)  # convierte a str el tipo de particula (endoB)
			lista_exoB.append(v)  # agrega valores de index de las particulas endoB a lista_endoB
			# lista_A_B_exoB[a1] = lista_exoB		#agrega al diccionario lesta la key endoB y les asigna sus values = lista con idex de tipo endoB
			# print ("lista endoB es ", lista_endoB)
			count = count + 1
		elif a == 'endo':
			lista_endo.append(count)
			count = count + 1
		elif a == 'endoB':
			# print(v)
			a1 = str(a)  # convierte a str el tipo de particula (endoB)
			lista_endoB.append(v)  # agrega valores de index de las particulas endoB a lista_endoB
			lista_coreB_endoB[
				a1] = lista_endoB  # agrega al diccionario lesta la key endoB y les asigna sus values = lista con idex de tipo endoB
			# print ("lista endoB es ", lista_endoB)
			lista_endoB_index.append(count)
			count = count + 1
		else:
			count = count +1
	if len(lista_endoB) == 0 and len(lista_endo) == 0 and len(lista_exo) == 0 and len(
			lista_exoB) == 0:  # no hay enzima unida
		recipe.change_topology_type("TopologyB")
	if len(lista_exoB) > 0 and len(lista_core_coreB) == 0 and len(vertices) > 3:  # exo unida a un frag de 2 (hA-hB)
		if len(lista_exoB) == 1:
			cant = len(lista_exoB)
			print(vertices)
			print("Este es cant", cant)
			print(lista_exoB, lista_A_B, lista_headA, lista_headB)
			reac = random.sample(range(0, len(lista_A_B)), cant)  # que par coreB-endoB van a reaccionar
			print(reac)
			for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
				print("i : " + str(i))
				z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
				print("z : " + str(z))
				for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
					if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
						if topology.particle_type_of_vertex(n) == "A":
							recipe.remove_edge(z, vertices.index(n))
							recipe.remove_edge(vertices.index(n), lista_headB[0])
							recipe.change_particle_type(z, "exo")
							recipe.change_particle_type(vertices.index(n), "celobiosa")
							recipe.change_particle_type(lista_headB[0], "celobiosaX")
						elif topology.particle_type_of_vertex(n) == "B":
							recipe.remove_edge(z, vertices.index(n))
							recipe.remove_edge(vertices.index(n), lista_headA[0])
							recipe.change_particle_type(z, "exo")
							recipe.change_particle_type(vertices.index(n), "celobiosa")
							recipe.change_particle_type(lista_headA[0], "celobiosaX")
		if len(lista_exoB) == 2:
			z = vertices.index(lista_exoB[0])  # index del exoB en la posicion 0 --> a reaccionar
			y = vertices.index(lista_exoB[1])  # index del exoB en la posicion 1 --> a liberarse
			print("z : " + str(z))
			print("y : " + str(y))
			for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
				if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
					a = n
					if topology.particle_type_of_vertex(a) == "A":
						for m in lista_A_B:
							if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(m):
								b = m
								recipe.remove_edge(z, vertices.index(a))
								recipe.remove_edge(y, vertices.index(b))
								recipe.remove_edge(vertices.index(a), vertices.index(b))
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(y, "exo")
								recipe.change_particle_type(vertices.index(a), "celobiosa")
								recipe.change_particle_type(vertices.index(b), "celobiosaX")
					if topology.particle_type_of_vertex(a) == "B":
						for m in lista_A_B:
							if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(m):
								b = m
								recipe.remove_edge(z, vertices.index(a))
								recipe.remove_edge(y, vertices.index(b))
								recipe.remove_edge(vertices.index(a), vertices.index(b))
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(y, "exo")
								recipe.change_particle_type(vertices.index(a), "celobiosa")
								recipe.change_particle_type(vertices.index(b), "celobiosaX")
	if len(lista_exoB) > 0 and len(lista_core_coreB) == 1:  # exo unida a un frag de 3 (hA-core-hB)
		print(vertices)
		if len(lista_A_B) == 1:  # una sola partic exo
			if len(lista_coreB_index) == 0:  # no hay particulas de endo unidas
				cant = 1
				print("Este es cant", cant)
				reac = random.sample(range(0, len(lista_A_B)), cant)  # que par coreB-endoB van a reaccionar
				print(reac)
				for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
					print("i : " + str(i))
					z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
					print("z : " + str(z))
					for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reacciona
						if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
							if topology.particle_type_of_vertex(n) == "A":
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_core[0])
								recipe.remove_edge(lista_core[0], lista_headB[0])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_core[0], "celobiosaX")
								recipe.change_particle_type(lista_headB[0], "glucosa")
							if topology.particle_type_of_vertex(n) == "B":
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_core[0])
								recipe.remove_edge(lista_core[0], lista_headA[0])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_core[0], "celobiosaX")
								recipe.change_particle_type(lista_headA[0], "glucosa")
			elif len(lista_coreB_index) == 1:  # hay particulas de endo unidas
				cant = 1
				print("Este es cant", cant)
				reac = random.sample(range(0, len(lista_A_B)), cant)  # que par coreB-endoB van a reaccionar
				for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
					print("i : " + str(i))
					z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
					print("z : " + str(z))
					for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
						if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
							if topology.particle_type_of_vertex(n) == "A":
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_coreB_index[0])
								recipe.remove_edge(lista_coreB_index[0], lista_headB[0])
								recipe.remove_edge(lista_coreB_index[0], lista_endoB_index[0])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_coreB_index[0], "celobiosaX")
								recipe.change_particle_type(lista_headB[0], "glucosa")
								recipe.change_particle_type(lista_endoB_index[0], "endo")
							if topology.particle_type_of_vertex(n) == "B":
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_coreB_index[0])
								recipe.remove_edge(lista_coreB_index[0], lista_headA[0])
								recipe.remove_edge(lista_coreB_index[0], lista_endoB_index[0])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_coreB_index[0], "celobiosaX")
								recipe.change_particle_type(lista_headA[0], "glucosa")
								recipe.change_particle_type(lista_endoB_index[0], "endo")
		elif len(lista_A_B) == 2:
			if len(lista_coreB_index) == 0:  # no hay particulas de endo unidas
				print("largo cel = 3, exoB = 2, endoB = 0")
				z = vertices.index(lista_exoB[0])  # index del exoB en la posicion 0 --> a reaccionar
				y = vertices.index(lista_exoB[1])  # index del exoB en la posicion 1 --> a liberarse
				print("z : " + str(z))
				print("y : " + str(y))
				for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
					if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
						a = n
						if topology.particle_type_of_vertex(a) == "A":
							for m in lista_A_B:
								if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(
										m):
									b = m
									recipe.remove_edge(vertices.index(b), y)
									recipe.remove_edge(vertices.index(a), z)
									recipe.remove_edge(lista_core[0], vertices.index(b))
									recipe.remove_edge(lista_core[0], vertices.index(a))
									recipe.change_particle_type(z, "exo")
									recipe.change_particle_type(y, "exo")
									recipe.change_particle_type(vertices.index(a), "celobiosa")
									recipe.change_particle_type(lista_core[0], "celobiosaX")
									recipe.change_particle_type(vertices.index(b), "glucosa")
						if topology.particle_type_of_vertex(n) == "B":
							for m in lista_A_B:
								if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(
										m):
									b = m
									recipe.remove_edge(vertices.index(b), y)
									recipe.remove_edge(vertices.index(a), z)
									recipe.remove_edge(lista_core[0], vertices.index(b))
									recipe.remove_edge(lista_core[0], vertices.index(a))
									recipe.change_particle_type(z, "exo")
									recipe.change_particle_type(y, "exo")
									recipe.change_particle_type(vertices.index(a), "celobiosa")
									recipe.change_particle_type(lista_core[0], "celobiosaX")
									recipe.change_particle_type(vertices.index(b), "glucosa")
			elif len(lista_coreB_index) == 1:  # hay particulas de endo unidas
				print("largo cel = 3, exoB = 2, endoB = 1")
				z = vertices.index(lista_exoB[0])  # index del exoB en la posicion 0 --> a reaccionar
				y = vertices.index(lista_exoB[1])  # index del exoB en la posicion 1 --> a liberarse
				print("z : " + str(z))
				print("y : " + str(y))
				for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
					if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
						a = n
						if topology.particle_type_of_vertex(a) == "A":
							for m in lista_A_B:
								if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(
										m):
									b = m
									recipe.remove_edge(vertices.index(b), y)
									recipe.remove_edge(vertices.index(a), z)
									recipe.remove_edge(lista_coreB_index[0], vertices.index(b))
									recipe.remove_edge(lista_coreB_index[0], vertices.index(a))
									recipe.remove_edge(lista_coreB_index[0], lista_endoB_index[0])
									recipe.change_particle_type(z, "exo")
									recipe.change_particle_type(y, "exo")
									recipe.change_particle_type(vertices.index(a), "celobiosa")
									recipe.change_particle_type(lista_coreB_index[0], "celobiosaX")
									recipe.change_particle_type(vertices.index(b), "glucosa")
									recipe.change_particle_type(lista_endoB_index[0], "endo")
						elif topology.particle_type_of_vertex(a) == "B":
							for m in lista_A_B:
								if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(
										m):
									b = m
									recipe.remove_edge(vertices.index(b), y)
									recipe.remove_edge(vertices.index(a), z)
									recipe.remove_edge(lista_coreB_index[0], vertices.index(b))
									recipe.remove_edge(lista_coreB_index[0], vertices.index(a))
									recipe.remove_edge(lista_coreB_index[0], lista_endoB_index[0])
									recipe.change_particle_type(z, "exo")
									recipe.change_particle_type(y, "exo")
									recipe.change_particle_type(vertices.index(a), "celobiosa")
									recipe.change_particle_type(lista_coreB_index[0], "celobiosaX")
									recipe.change_particle_type(vertices.index(b), "glucosa")
									recipe.change_particle_type(lista_endoB_index[0], "endo")

	if len(lista_exoB) > 0 and len(lista_core_coreB) == 2:  # exo unida a un frag de 4 (hA-core-core-hB)
		print(vertices)
		if len(lista_A_B) == 1:
			if len(lista_coreB_index) == 0:
				cant = 1
				print("Este es cant", cant)
				reac = random.sample(range(0, len(lista_A_B)),
									 cant)  # que par coreB-endoB van a reaccionar (tantos valores como determine variable cant). es un n que va de 0 a el valor de la cant de coreB menos 1 (menos 1 por que corrijo por el 0)
				print(reac)
				for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
					print("i : " + str(i))
					z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
					print("z : " + str(z))
					for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
						if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
							if topology.particle_type_of_vertex(n) == "A":
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_core[0])
								recipe.remove_edge(lista_core[0], lista_core[1])
								recipe.remove_edge(lista_core[1], lista_headB[0])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_core[0], "celobiosaX")
								recipe.change_particle_type(lista_core[1], "celobiosaX")
								recipe.change_particle_type(lista_headB[0], "celobiosa")
							elif topology.particle_type_of_vertex(n) == "B":
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_core[1])
								recipe.remove_edge(lista_core[1], lista_core[0])
								recipe.remove_edge(lista_core[0], lista_headA[0])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_core[1], "celobiosaX")
								recipe.change_particle_type(lista_core[0], "celobiosaX")
								recipe.change_particle_type(lista_headA[0], "celobiosa")
			elif len(lista_coreB_index) > 0:
				print("Atasco!")
		elif len(lista_A_B) == 2:
			if len(lista_coreB_index) == 0:
				print("largo cel = 4, exoB = 2")
				z = vertices.index(lista_exoB[0])  # index del exoB en la posicion 0 --> a reaccionar
				y = vertices.index(lista_exoB[1])  # index del exoB en la posicion 1 --> a liberarse
				print("z : " + str(z))
				print("y : " + str(y))
				for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
					if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
						a = n
						if topology.particle_type_of_vertex(a) == "A":
							for m in lista_A_B:
								if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(
										m):
									b = m
									recipe.remove_edge(vertices.index(b), y)
									recipe.remove_edge(vertices.index(a), z)
									recipe.remove_edge(vertices.index(a), lista_core[0])
									recipe.remove_edge(lista_core[0], lista_core[1])
									recipe.remove_edge(lista_core[1], vertices.index(b))
									recipe.change_particle_type(z, "exo")
									recipe.change_particle_type(y, "exo")
									recipe.change_particle_type(vertices.index(a), "celobiosa")
									recipe.change_particle_type(lista_core[0], "celobiosaX")
									recipe.change_particle_type(lista_core[1], "celobiosaX")
									recipe.change_particle_type(vertices.index(b), "celobiosa")
						elif topology.particle_type_of_vertex(a) == "B":
							for m in lista_A_B:
								if str(vertices[y].neighbors()[0]) == str(m) or str(vertices[y].neighbors()[1]) == str(
										m):
									b = m
									recipe.remove_edge(vertices.index(b), y)
									recipe.remove_edge(vertices.index(a), z)
									recipe.remove_edge(vertices.index(b), lista_core[0])
									recipe.remove_edge(lista_core[0], lista_core[1])
									recipe.remove_edge(lista_core[1], vertices.index(a))
									recipe.change_particle_type(z, "exo")
									recipe.change_particle_type(y, "exo")
									recipe.change_particle_type(vertices.index(a), "celobiosa")
									recipe.change_particle_type(lista_core[0], "celobiosaX")
									recipe.change_particle_type(lista_core[1], "celobiosaX")
									recipe.change_particle_type(vertices.index(b), "celobiosa")
			elif len(lista_coreB_index) > 0:
				print("Atasco 1")
	if len(lista_exoB) > 0 and len(lista_core_coreB) > 2:  # exo unida a un frag de > de 4 (hA-(core)n-hB) con n > 2
		print(vertices)
		if len(lista_coreB_index) == 0:
			# cant = random.randint(1, len(lista_A_B))
			cant = 1
			print("Este es cant", cant)
			reac = random.sample(range(0, len(lista_A_B)),
								 cant)  # que par coreB-endoB van a reaccionar (tantos valores como determine variable cant). es un n que va de 0 a el valor de la cant de coreB menos 1 (menos 1 por que corrijo por el 0)
			print(reac)
			for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
				print("i : " + str(i))
				z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
				print("z : " + str(z))
				for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
					if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
						print("n : " + str(n))
						print(topology.particle_type_of_vertex(n))
						if topology.particle_type_of_vertex(n) == "A":
							print("VER1: " + str(vertices.index(n)) + str(lista_core[0]) + str(lista_core[1]))
							recipe.remove_edge(vertices.index(n), z)
							recipe.remove_edge(vertices.index(n), lista_core[0])
							recipe.remove_edge(lista_core[0], lista_core[1])
							recipe.change_particle_type(z, "exo")
							recipe.change_particle_type(vertices.index(n), "celobiosa")
							recipe.change_particle_type(lista_core[0], "celobiosaX")
							recipe.change_particle_type(lista_core[1], "headA")
						elif topology.particle_type_of_vertex(n) == "B":
							print("VER2: " + str(vertices.index(n)) + str(lista_core[0]) + str(lista_core[1]))
							recipe.remove_edge(vertices.index(n), z)
							recipe.remove_edge(vertices.index(n), lista_core[-1])
							recipe.remove_edge(lista_core[-1], lista_core[-2])
							recipe.change_particle_type(z, "exo")
							recipe.change_particle_type(vertices.index(n), "celobiosa")
							recipe.change_particle_type(lista_core[-1], "celobiosaX")
							recipe.change_particle_type(lista_core[-2], "headB")
		elif len(lista_coreB_index) > 0:
			# cant = random.randint(1, len(lista_A_B))
			cant = 1
			print("Este es cant", cant)
			reac = random.sample(range(0, len(lista_A_B)),
								 cant)  # que par coreB-endoB van a reaccionar (tantos valores como determine variable cant). es un n que va de 0 a el valor de la cant de coreB menos 1 (menos 1 por que corrijo por el 0)
			print(reac)
			for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
				print("i : " + str(i))
				z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
				print("z : " + str(z))
				for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
					if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
						print("n : " + str(n))
						print(topology.particle_type_of_vertex(n))
						if topology.particle_type_of_vertex(n) == "A":
							condA = lista_core_coreB[0] in lista_coreB_index or lista_core_coreB[1] in lista_coreB_index
							if condA == False:
								print("VER3: " + str(vertices.index(n)) + str(lista_core[0]) + str(lista_core[1]))
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_core[0])
								recipe.remove_edge(lista_core[0], lista_core[1])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_core[0], "celobiosaX")
								recipe.change_particle_type(lista_core[1], "headA")
							elif condA == True:
								print("Atasco 2")
						elif topology.particle_type_of_vertex(n) == "B":
							condB = lista_core_coreB[-1] in lista_coreB_index or lista_core_coreB[
								-2] in lista_coreB_index
							if condB == False:
								print("VER4: " + str(vertices.index(n)) + str(lista_core[-1]) + str(lista_core[-2]))
								recipe.remove_edge(vertices.index(n), z)
								recipe.remove_edge(vertices.index(n), lista_core[-1])
								recipe.remove_edge(lista_core[-1], lista_core[-2])
								recipe.change_particle_type(z, "exo")
								recipe.change_particle_type(vertices.index(n), "celobiosa")
								recipe.change_particle_type(lista_core[-1], "celobiosaX")
								recipe.change_particle_type(lista_core[-2], "headB")
							elif condB == True:
								print("Atasco 3")
	return recipe


system.topologies.add_structural_reaction(name="exog1",
										  topology_type="TopologyC",
										  reaction_function=exog1_reaction_function,
										  rate_function=exog1_rate_function)


def exog2_rate_function(topology):  # k-1 disociacion exoglucanasa
	return 0.01


def exog2_reaction_function(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertices = topology.get_graph().get_vertices()
	edges = topology.get_graph().get_edges()
	lista_core = []  # lista con index de core en orden
	lista_coreB_index = []  # lista con index de coreB en orden
	lista_coreB = []
	lista_endoB = []  # lista con index de endoB en orden
	lista_endo = []
	lista_exoB = []  # lista con index de endoB en orden
	lista_exo = []
	lista_core_coreB = []  # lista con index de core y coreB en orden
	lista_coreB_endoB = {}  # diccionario con keys = nombre de particula, values = listas con index de particulas en orden
	lista_headB = []  # lista con index de headB en orden
	lista_headA = []
	lista_A_B = []
	lista_A_B_exoB = []
	count = 0
	for v in vertices:
		a = topology.particle_type_of_vertex(v)
		print(a, count, topology.particle_id_of_vertex(v))
		if a == 'headB':
			lista_headB.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headB es", lista_headB)
			largocadena = count + 1  # calculo de largo de cadena. Sirve para limitar la reaccion a cierto largo de cadena
			anteultimo = largocadena - 2  # calculo de anteultimo valor de cadena. Sirve para que endo no reaccione con borde generando glucosa.
			# print("largo cadena es " + str(largocadena))
			count = count + 1
		elif a == 'headA':
			lista_headA.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headA es", lista_headA)
			count = count + 1
		elif a == 'A':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_A_B.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			# lista_A_index.append(count)
			# lista_A_B_exoB[a0] = lista_A_B		#agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'B':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_A_B.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			# lista_A_index.append(count)
			# lista_core_coreB.append(count)
			# lista_A_B_exoB[a0] = lista_A_B		#agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'core':
			lista_core.append(count)  # agrega valores de index de las particulas core a lista_core
			lista_core_coreB.append(count)
			# print("lista core es ", lista_core)
			count = count + 1
		elif a == 'coreB':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_coreB.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			lista_coreB_index.append(count)
			lista_core_coreB.append(count)
			lista_coreB_endoB[
				a0] = lista_coreB  # agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
	#	elif a == 'F':
	#		count = count + 1
	#	elif a == 'F1':
	#		count = count + 1
	#	elif a == 'F2':
	#		count = count + 1
	#	elif a == 'F3':
	#		count = count + 1
	#	elif a == 'F4':
	#		count = count + 1
	#	elif a == 'F5':
	#		count = count + 1
	#	elif a == 'F6':
	#		count = count + 1
	#	elif a == 'F7':
	#		count = count + 1
	#	elif a == 'F8':
	#		count = count + 1
	#	elif a == 'F9':
	#		count = count + 1
	#	elif a == 'F10':
	#		count = count + 1
	#	elif a == 'F11':
	#		count = count + 1
	#	elif a == 'F12':
	#		count = count + 1
	#	elif a == 'F13':
	#		count = count + 1
	#	elif a == 'F14':
	#		count = count + 1
	#	elif a == 'F15':
	#		count = count + 1
	#	elif a == 'F16':
	#		count = count + 1
	#	elif a == 'F17':
	#		count = count + 1
	#	elif a == 'F18':
	#		count = count + 1
	#	elif a == 'F19':
	#		count = count + 1
	#	elif a == 'F20':
	#		count = count + 1
		elif a == 'exo':
			lista_exo.append(count)
			count = count + 1
		elif a == 'exoB':
			# print(v)
			a1 = str(a)  # convierte a str el tipo de particula (endoB)
			lista_exoB.append(v)  # agrega valores de index de las particulas endoB a lista_endoB
			# lista_A_B_exoB[a1] = lista_exoB		#agrega al diccionario lesta la key endoB y les asigna sus values = lista con idex de tipo endoB
			# print ("lista endoB es ", lista_endoB)
			count = count + 1
		elif a == 'endo':
			lista_endo.append(count)
			count = count + 1
		elif a == 'endoB':
			# print(v)
			a1 = str(a)  # convierte a str el tipo de particula (endoB)
			lista_endoB.append(v)  # agrega valores de index de las particulas endoB a lista_endoB
			lista_coreB_endoB[
				a1] = lista_endoB  # agrega al diccionario lesta la key endoB y les asigna sus values = lista con idex de tipo endoB
			# print ("lista endoB es ", lista_endoB)
			count = count + 1
		else:
			count = count +1
	if len(lista_endoB) == 0 and len(lista_endo) == 0 and len(lista_exo) == 0 and len(
			lista_exoB) == 0:  # no hay enzima unida
		recipe.change_topology_type("TopologyB")
	if len(lista_exoB) > 0:
		print(vertices)
		cant = 1
		print("Este es cant", cant)
		reac = random.sample(range(0, len(lista_A_B)),
							 cant)  # que par coreB-endoB van a reaccionar (tantos valores como determine variable cant). es un n que va de 0 a el valor de la cant de coreB menos 1 (menos 1 por que corrijo por el 0)
		print(reac)
		for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
			print("i : " + str(i))
			z = vertices.index(lista_exoB[i])  # index del exoB en la posicion i de reac.
			print("z : " + str(z))
			for n in lista_A_B:  # ubicar cual de los vecinos de exoB es la particula a reaccionar
				if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
					if topology.particle_type_of_vertex(n) == "A":
						recipe.remove_edge(n, lista_exoB[i])
						recipe.change_particle_type(n, "headA")
						recipe.change_particle_type(lista_exoB[i], "exo")
					elif topology.particle_type_of_vertex(n) == "B":
						recipe.remove_edge(n, lista_exoB[i])
						recipe.change_particle_type(n, "headB")
						recipe.change_particle_type(lista_exoB[i], "exo")
	return recipe


system.topologies.add_structural_reaction(name="exog2",
										  topology_type="TopologyC",
										  reaction_function=exog2_reaction_function,
										  rate_function=exog2_rate_function)

## endoglucanase ##


system.topologies.add_spatial_reaction(
	'endo bindA: TopologyB(core) + TopologyC(endo) -> TopologyC(coreB--endoB) [self=true]', rate=0.1, radius=endo + core
	# 100 sM-1 (k1) binding of endo to core (cellulose topology B)
)

system.topologies.add_spatial_reaction(
	'endo bindAA: TopologyC(core) + TopologyC(endo) -> TopologyC(coreB--endoB) [self=true]', rate=0.1, radius=endo + core
	# (k1) binding of endo to core (cellulose topology C)
)  # 150


def endog1_rate_function(topology):  # kcat para reaccion de endoglucanasa
	return 0.001


def endog1_reaction_function(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertices = topology.get_graph().get_vertices()
	edges = topology.get_graph().get_edges()
	lista_core = []  # lista con index de core en orden
	lista_coreB_index = []  # lista con index de coreB en orden
	lista_coreB = []
	lista_endoB = []  # lista con index de endoB en orden
	lista_endo = []
	lista_exo = []
	lista_exoB = []
	lista_core_coreB = []  # lista con index de core y coreB en orden
	lista_coreB_endoB = {}  # diccionario con keys = nombre de particula, values = listas con index de particulas en orden
	lista_headB = []  # lista con index de headB en orden
	lista_headA = []
	lista_A_B = []
	count = 0
	for v in vertices:
		# print ("vertices ", vertices[0])					!!!!!!!
		# v1 = v1_ref.get()
		# v2 = v2_ref.get()
		# c = v.neighbors()[-1].get().particle_index
		# print(v.neighbors()[0])
		# print (v1.neighbors())
		a = topology.particle_type_of_vertex(v)
		print(a, count, topology.particle_id_of_vertex(v))
		# print (v)
		# vec1 = str(vertices[1].neighbors()[0])
		# print (vec1)
		# print (str(v))
		# print (str(v) == vec1)
		if a == 'headB':
			lista_headB.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headB es", lista_headB)
			largocadena = count + 1  # calculo de largo de cadena. Sirve para limitar la reaccion a cierto largo de cadena
			anteultimo = largocadena - 2  # calculo de anteultimo valor de cadena. Sirve para que endo no reaccione con borde generando glucosa.
			# print("largo cadena es " + str(largocadena))
			count = count + 1
		elif a == 'headA':
			lista_headA.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headA es", lista_headA)
			count = count + 1
		elif a == 'core':
			lista_core.append(count)  # agrega valores de index de las particulas core a lista_core
			lista_core_coreB.append(count)
			# print("lista core es ", lista_core)
			count = count + 1
		elif a == 'coreB':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_coreB.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			lista_coreB_index.append(count)
			lista_core_coreB.append(count)
			lista_coreB_endoB[
				a0] = lista_coreB  # agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'exo':
			lista_exo.append(count)
			count = count + 2
		elif a == 'exoB':
			lista_exoB.append(count)
			count = count + 2
		elif a == 'A':
			lista_A_B.append(count)
			count = count + 1
		elif a == 'B':
			lista_A_B.append(count)
			count = count + 1
		elif a == 'endo':
			lista_endo.append(count)
			count = count + 1
		elif a == 'endoB':
			# print(v)
			a1 = str(a)  # convierte a str el tipo de particula (endoB)
			lista_endoB.append(v)  # agrega valores de index de las particulas endoB a lista_endoB
			lista_coreB_endoB[
				a1] = lista_endoB  # agrega al diccionario lesta la key endoB y les asigna sus values = lista con idex de tipo endoB
			# print ("lista endoB es ", lista_endoB)
			count = count + 1
		else:
			count = count + 1
	if len(lista_endoB) == 0 and len(lista_endo) == 0 and len(lista_exo) == 0 and len(lista_exoB) == 0:
		recipe.change_topology_type("TopologyB")
	if len(lista_coreB) > 0:
		print(vertices)
		cant = random.randint(1, len(lista_coreB))
		print("Este es cant", cant)
		reac = random.sample(range(0, len(lista_coreB)),
							 cant)  # que par coreB-endoB van a reaccionar (tantos valores como determine variable cant). es un n que va de 0 a el valor de la cant de coreB menos 1 (menos 1 por que corrijo por el 0)
		print(reac)
		lado = random.randint(0, 1)
		for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
			print(i)
			z = vertices.index(lista_coreB_endoB["endoB"][i])
			print(z)
			for n in lista_coreB_endoB["coreB"]:  # ubicar cual es el coreB unido a endoB
				print(n)
				print(vertices[z].neighbors()[0], vertices[z].neighbors()[1])
				if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
					cond1 = lista_core_coreB.index(vertices.index(n)) >= 1 and lista_core_coreB.index(
						vertices.index(n)) <= len(
						lista_core_coreB) - 2  # cond1 chequea que el coreB a reaccionar este separado de al menos 1 particula antes de los head (>= 1 por que el 0 de esta lista es el primer core que hay y -2 por que len lista da un numero que no arranca a contar de 0 por lo tanto hay que restarle 1 y 1 menos ya que la ultima particula de la lista es el core anterior a headB).
					if cond1 is True:
						cond2 = lista_core_coreB[lista_core_coreB.index(vertices.index(n)) + 1] in lista_coreB_index or \
								lista_core_coreB[lista_core_coreB.index(vertices.index(
									n)) - 1] in lista_coreB_index  # no tiene otro coreB pegado en +1 o -1
						if cond2 is False and lado == 1:
							print("reaccion1.1", n, lista_coreB_endoB["endoB"][i], vertices.index(n),
								  lista_core_coreB[lista_core_coreB.index(vertices.index(n)) + 1])
							recipe.remove_edge(n, lista_coreB_endoB["endoB"][
								i])  # corte enlace entre indice coreB y indice endoB (separa enzima)
							recipe.remove_edge(vertices.index(n), lista_core_coreB[lista_core_coreB.index(
								vertices.index(
									n)) + 1])  # corte enlace entre indice coreB y indice coreB +1 (separa celulosa)
							recipe.change_particle_type(lista_core_coreB[lista_core_coreB.index(vertices.index(n)) + 1],
														"headA")  # conversion particula indice coreB+1 en extremo headA
							recipe.change_particle_type(n,
														"headB")  # conversion particula indice coreB en extremo headB
							recipe.change_particle_type(lista_coreB_endoB["endoB"][i],
														"endo")  # conversion particula indice endoB en endo (regeneracion enzima)
						elif cond2 is False and lado == 0:  # si la particula que le sigue en index es core, entonces el enlace a romper es entre index de coreB e index coreB +1
							print("reaccion1.0", n, lista_coreB_endoB["endoB"][i], vertices.index(n),
								  lista_core_coreB[lista_core_coreB.index(vertices.index(n)) - 1])
							recipe.remove_edge(n, lista_coreB_endoB["endoB"][
								i])  # corte enlace entre indice coreB y indice endoB (separa enzima)
							recipe.remove_edge(vertices.index(n), lista_core_coreB[lista_core_coreB.index(
								vertices.index(
									n)) - 1])  # corte enlace entre indice coreB y indice coreB +1 (separa celulosa)
							recipe.change_particle_type(lista_core_coreB[lista_core_coreB.index(vertices.index(n)) - 1],
														"headB")  # conversion particula indice coreB+1 en extremo headA
							recipe.change_particle_type(n,
														"headA")  # conversion particula indice coreB en extremo headB
							recipe.change_particle_type(lista_coreB_endoB["endoB"][i],
														"endo")  # conversion particula indice endoB en endo (regeneracion enzima)
						elif cond2 is True:
							print("No hay reaccion")
							recipe.remove_edge(n, lista_coreB_endoB["endoB"][
								i])  # corte enlace entre indice coreB y indice endoB (separa enzima)
							recipe.change_particle_type(n, "core")  # conversion particula indice coreB en extremo core
							recipe.change_particle_type(lista_coreB_endoB["endoB"][i],
														"endo")  # conversion particula indice endoB en endo (regeneracion enzima)
					elif cond1 is False:
						print("No hay reaccion")
						recipe.remove_edge(n, lista_coreB_endoB["endoB"][
							i])  # corte enlace entre indice coreB y indice endoB (separa enzima)
						recipe.change_particle_type(n, "core")  # conversion particula indice coreB en extremo core
						recipe.change_particle_type(lista_coreB_endoB["endoB"][i],
													"endo")  # conversion particula indice endoB en endo (regeneracion enzima)
	return recipe


system.topologies.add_structural_reaction(name="endog1",
										  topology_type="TopologyC",
										  reaction_function=endog1_reaction_function,
										  rate_function=endog1_rate_function)


def endog2_rate_function(topology):  # koff para reaccion de endoglucanasa
	return 0.002


def endog2_reaction_function(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertices = topology.get_graph().get_vertices()
	edges = topology.get_graph().get_edges()
	lista_core = []  # lista con index de core en orden
	lista_coreB_index = []  # lista con index de coreB en orden
	lista_coreB = []
	lista_endo = []
	lista_endoB = []  # lista con index de endoB en orden
	lista_exo = []
	lista_exoB = []
	lista_core_coreB = []  # lista con index de core y coreB en orden
	lista_coreB_endoB = {}  # diccionario con keys = nombre de particula, values = listas con index de particulas en orden
	lista_headB = []  # lista con index de headB en orden
	lista_headA = []
	lista_A_B = []
	count = 0
	for v in vertices:
		a = topology.particle_type_of_vertex(v)
		print(a, count, topology.particle_id_of_vertex(v))
		if a == 'headB':
			lista_headB.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headB es", lista_headB)
			largocadena = count + 1  # calculo de largo de cadena. Sirve para limitar la reaccion a cierto largo de cadena
			anteultimo = largocadena - 2  # calculo de anteultimo valor de cadena. Sirve para que endo no reaccione con borde generando glucosa.
			# print("largo cadena es " + str(largocadena))
			count = count + 1
		elif a == 'headA':
			lista_headA.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headA es", lista_headA)
			count = count + 1
		elif a == 'A':
			lista_A_B.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headA es", lista_headA)
			count = count + 1
		elif a == 'B':
			lista_A_B.append(count)  # agrega valores de index de las particulas headB a lista_headB
			# print("lista headA es", lista_headA)
			count = count + 1
		elif a == 'core':
			lista_core.append(count)  # agrega valores de index de las particulas core a lista_core
			lista_core_coreB.append(count)
			# print("lista core es ", lista_core)
			count = count + 1
		elif a == 'coreB':
			# print(v)
			a0 = str(a)  # convierte a str el tipo de particula (coreB)
			lista_coreB.append(v)  # agrega valores de index de las particulas coreB a lista_coreB
			lista_coreB_index.append(count)
			lista_core_coreB.append(count)
			lista_coreB_endoB[
				a0] = lista_coreB  # agrega al diccionario lesta la key coreB y les asigna sus values = lista con idex de tipo coreB
			# print ("lista coreB es ", lista_coreB)
			count = count + 1
		elif a == 'endo':
			lista_endo.append(count)
			count = count + 1
		elif a == 'exo':
			lista_exo.append(count)
			count = count + 1
		elif a == 'exoB':
			lista_exoB.append(count)
			count = count + 1
		elif a == 'endoB':
			# print(v)
			a1 = str(a)  # convierte a str el tipo de particula (endoB)
			lista_endoB.append(v)  # agrega valores de index de las particulas endoB a lista_endoB
			lista_coreB_endoB[
				a1] = lista_endoB  # agrega al diccionario lesta la key endoB y les asigna sus values = lista con idex de tipo endoB
			# print ("lista endoB es ", lista_endoB)
			count = count + 1
		else:
			count = count + 1
	if len(lista_endoB) == 0 and len(lista_endo) == 0 and len(lista_exo) == 0 and len(lista_exoB) == 0:
		recipe.change_topology_type("TopologyB")
	if len(lista_endoB) > 0:
		print(vertices)
		cant = random.randint(1,
							  len(lista_coreB))  # cant = seleccion al azar de cantidad de reacciones a ocurrir (n entre 1 y cant de coreB)
		print("este es cant")
		print(cant)
		reac = random.sample(range(0, len(lista_coreB)),
							 cant)  # que par coreB-endoB van a reaccionar (tantos valores como determine variable cant, AHORA YA NO ES ASI POR QUE DABA ERROR CON LA ACTUALIZACION DE CODIGO). es un n que va de 0 a el valor de la cant de coreB menos 1 (menos 1 por que corrijo por el 0)
		# print("este es reac")
		for i in reac:  # los i van a ser iguales a los valores seleccionados en reac y como fue aleatorio no estaran en orden
			z = vertices.index(lista_coreB_endoB["endoB"][i])
			for n in lista_coreB_endoB["coreB"]:
				# print (vertices.index(lista_coreB_endoB["coreB"][i]))
				if str(vertices[z].neighbors()[0]) == str(n) or str(vertices[z].neighbors()[1]) == str(n):
					# print("TRUE")
					recipe.remove_edge(n, lista_coreB_endoB["endoB"][
						i])  # corte enlace entre indice coreB y indice endoB (separa enzima)
					recipe.change_particle_type(n, "core")  # conversion particula indice coreB en extremo headB
					recipe.change_particle_type(lista_coreB_endoB["endoB"][i],
												"endo")  # conversion particula indice endoB en endo (regeneracion enzima)
					print("disociacion")
	return recipe


system.topologies.add_structural_reaction(name="endog2",
										  topology_type="TopologyC",
										  reaction_function=endog2_reaction_function,
										  rate_function=endog2_rate_function)


## remove leftover cellulose ##

def decayrestos_rate_function(topology):
	return 0.3


def decayrestos_reaction_function(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertices = topology.get_graph().get_vertices()
	edges = topology.get_graph().get_edges()
	nvertices = len(vertices)
	lista_core = []  # lista con index de core en orden
	lista_coreB = []  # lista con index de coreB en orden
	lista_endoB = []  # lista con index de endoB en orden
	lista_coreB_endoB = {}  # diccionario con keys = nombre de particula, values = listas con index de particulas en orden
	lista_headB = []  # lista con index de headB en orden
	for v in vertices:
		b = topology.particle_type_of_vertex(v)
		# print(a,b)
		if len(vertices) == 1 and b == 'headA':
			print("restos")
			recipe.change_particle_type(v, "glucosa")
		if len(vertices) == 1 and b == 'headB':
			print("restos")
			recipe.change_particle_type(v, "glucosa")
	for e in edges:
		v1_ref, v2_ref = e[0], e[1]
		v1 = v1_ref.get()
		v2 = v2_ref.get()
		a = [topology.particle_type_of_vertex(v1), topology.particle_type_of_vertex(v2)]
		b = [v1.particle_index, v2.particle_index]
		if len(vertices) == 2:
			if a[0] == 'headA' and a[1] == 'headA':
				print("restos")
				recipe.separate_vertex(0)
				recipe.change_particle_type(0, "celobiosa")
				recipe.change_particle_type(1, "celobiosaX")
			elif a[0] == 'headA' and a[1] == 'headB':
				print("restos")
				recipe.separate_vertex(0)
				recipe.change_particle_type(0, "celobiosa")
				recipe.change_particle_type(1, "celobiosaX")
			elif a[0] == 'headB' and a[1] == 'headA':
				print("restos")
				recipe.separate_vertex(0)
				recipe.change_particle_type(0, "celobiosa")
				recipe.change_particle_type(1, "celobiosaX")
			elif a[0] == 'headB' and a[1] == 'headB':
				print("restos")
				recipe.separate_vertex(0)
				recipe.change_particle_type(0, "celobiosa")
				recipe.change_particle_type(1, "celobiosaX")
	return recipe


system.topologies.add_structural_reaction(name="restos",
										  topology_type="TopologyB",
										  reaction_function=decayrestos_reaction_function,
										  rate_function=decayrestos_rate_function)

### SIMULACION ###

simulation = system.simulation(kernel="CPU")
simulation.kernel_configuration.n_threads = 1
simulation.output_file = 'output'
simulation.reaction_handler = "UncontrolledApproximation"
simulation.observe.particle_positions(stride=1000, types =["F1"])
simulation.observe.number_of_particles(stride=1000, types=["celobiosa", "celobiosaB", "celobiosaI", "glucosa", "headB", "B", "endo", "exo"], callback=lambda n: print("celobiosa, celobiosaB, celobiosaI, glucosa, headB, B:", n))



###Funciones de la simulacion###
init_pos = []
msd_pos = []
def simular_BLS(enzima, cantidad = n*10):
    """Agrega la cantidad ya definida de decameros de BLS"""
#    if enzima == "endo":
#        e1 = "endo"
#        e2 = "endo"
#    elif enzima == "exo":
#        e1 = "exo"
#        e2 = "exo"
#    elif enzima == "endo_exo":
#        e1 = "endo"
#        e2 = "exo"
    count = 0
    while count < cantidad:
        x = (50 * np.random.random(1)) - 25
        x = x.item()
        y = (50 * np.random.random(1)) - 25
        y = y.item()
        z = (50 * np.random.random(1)) - 25
        z = z.item()
        if (count / 10 % 2) == 0:
            e1 = "endo"
            e2 = "endo"
        else:
            e1 = "exo"
            e2 = "exo"
        if int(x) in range(-18, 18) and int(y) in range(-1, 14) and int(z) in range(-3, 12):
                count = count
        else:
                posiciones = [[0. + x, 0. + y, 0. + z], [2. + x, 0. + y, 0. + z], [2.61803 + x, 1.9021 + y, 0. + z], [1. + x, 3.07767 + y, 0. + z], [-0.61803 + x, 1.9021 + y, 0. + z], [0. + x, 0. + y, 2. + z], [2. + x, 0. + y, 2. + z], [2.61803 + x, 1.9021 + y, 2. + z], [1. + x, 3.07767 + y, 2. + z],[-0.61803 + x, 1.9021 + y, 2. + z]]
                init_pos.append(posiciones)
                msd_pos.append(posiciones[0])

                BLS = simulation.add_topology(topology_type="TopologyA", particle_types=BLS_list[count:count+10],positions=np.array(posiciones))
                graph = BLS.get_graph()
                # cierre de pentameros
                graph.add_edge(0, 1)
                graph.add_edge(1, 2)
                graph.add_edge(2, 3)
                graph.add_edge(3, 4)
                graph.add_edge(4, 0)
                graph.add_edge(5, 6)
                graph.add_edge(6, 7)
                graph.add_edge(7, 8)
                graph.add_edge(8, 9)
                graph.add_edge(9, 5)
                # union entre pentameros
                graph.add_edge(0, 5)
                graph.add_edge(1, 6)
                graph.add_edge(2, 7)
                graph.add_edge(3, 8)
                graph.add_edge(4, 9)
                print(BLS_list[count:count:+10])

                posiciones1 = [[0. + x, 0. + y, -1.3+ z], [0. + x, 0. + y, -2.6 + z]]
                enz1 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count], e1],positions=np.array(posiciones1))  # 1
                print(F_list[count])
                graph = enz1.get_graph()
                graph.add_edge(0, 1)

                posiciones2 = [[2. + x, 0. + y, -1.3 + z], [2. + x, 0. + y, -2.6 + z]]
                enz2 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+1], e2],positions=np.array(posiciones2))  # 2
                print(F_list[count+1])
                graph = enz2.get_graph()
                graph.add_edge(0, 1)

                posiciones3 = [[2.61803 + x, 1.9021 + y, -1.3 + z], [2.61803 + x, 1.9021 + y, -2.6 + z]]
                enz3 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+2], e1],positions=np.array(posiciones3))  # 3
                print(F_list[count + 2])
                graph = enz3.get_graph()
                graph.add_edge(0, 1)

                posiciones4 = [[1. + x, 3.07767 + y, -1.3 + z], [1. + x, 3.07767 + y, -2.6 + z]]
                enz4 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+3], e2],positions=np.array(posiciones4))  # 4
                print(F_list[count + 3])
                graph = enz4.get_graph()
                graph.add_edge(0, 1)

                posiciones5 = [[-0.61803 + x, 1.9021 + y, -1.3 + z], [-0.61803 + x, 1.9021 + y, -2.6 + z]]
                enz5 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+4], e1],positions=np.array(posiciones5))  # 5
                print(F_list[count + 4])
                graph = enz5.get_graph()
                graph.add_edge(0, 1)

                posiciones6 = [[0. + x, 0. + y, 3.3 + z], [0. + x, 0. + y, 4.6 + z]]
                enz6 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+5], e2],positions=np.array(posiciones6))  # 6
                print(F_list[count + 5])
                graph = enz6.get_graph()
                graph.add_edge(0, 1)

                posiciones7 = [[2. + x, 0. + y, 3.3 + z], [2. + x, 0. + y, 4.6 + z]]
                enz7 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+6], e1],positions=np.array(posiciones7))  # 7
                print(F_list[count + 6])
                graph = enz7.get_graph()
                graph.add_edge(0, 1)

                posiciones8 = [[2.61803 + x, 1.9021 + y, 3.3 + z], [2.61803 + x, 1.9021 + y, 4.6 + z]]
                enz8 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+7], e2],positions=np.array(posiciones8))  # 8
                print(F_list[count + 7])
                graph = enz8.get_graph()
                graph.add_edge(0, 1)

                posiciones9 = [[1. + x, 3.07767 + y, 3.3 + z], [1. + x, 3.07767 + y, 4.6 + z]]
                enz9 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+8], e1],positions=np.array(posiciones9))  # 9
                print(F_list[count + 8])
                graph = enz9.get_graph()
                graph.add_edge(0, 1)

                posiciones10 = [[-0.61803 + x, 1.9021 + y, 3.3 + z], [-0.61803 + x, 1.9021 + y, 4.6 + z]]
                enz10 = simulation.add_topology(topology_type="TopologyC", particle_types=[F_list[count+9], e2],positions=np.array(posiciones10))  # 10
                print(F_list[count + 9])
                graph = enz10.get_graph()
                graph.add_edge(0, 1)
                count += 10

def simular_LIBRES(enzima,cantidad=n*10):
	cantidad_endo = 0
	cantidad_exo = 0
	if enzima == "endo":
		cantidad_endo = cantidad
	elif enzima == "exo":
		cantidad_exo = cantidad
	elif enzima == "endo_exo":
		cantidad_exo = cantidad
		cantidad_endo = cantidad
	else:
		print("Tipo de enzima invalido. Utilizar: endo/exo/endo_exo")
	## ENDO ##
	count_endo = 0
	while count_endo < cantidad_endo:
		x = random.randint(-30, 30)
		y = random.randint(-30, 30)
		z = random.randint(-30, 30)

#		for position in range(len(ocupadas)):
#			if ocupadas[position][0] == x or ocupadas[position][1] == y or ocupadas[position][2] == z:
			#	count_endo = count_endo
			
		if int(x) in range(-18, 18) and int(y) in range(-1, 14) and int(z) in range(-3, 12):
			count_endo = count_endo
				
		else:
			count_endo += 1
			posiciones = [[0. + x, 0. + y, 0. + z], [1. + x, 0. + y, 0. + z]]
			msd_pos.append(posiciones[0])
			BLS = simulation.add_topology(topology_type="TopologyC", particle_types=["F1", "endo"],
										  positions=np.array(posiciones))
			graph = BLS.get_graph()
			graph.add_edge(0, 1)

	## EXO ##
	count_exo = 0
	while count_exo < cantidad_exo:
		x = random.randint(-30, 30)
		y = random.randint(-30, 30)
		z = random.randint(-30, 30)

		#for position in range(len(ocupadas)):
		#	if ocupadas[position][0] == x or ocupadas[position][1] == y or ocupadas[position][2] == z:
		#		count_exo = count_exo
			
		if int(x) in range(-18, 18) and int(y) in range(-1, 14) and int(z) in range(-3, 12):
			count_exo = count_exo
				
		else:
			count_exo += 1
			posiciones = [[0. + x, 0. + y, 0. + z], [1. + x, 0. + y, 0. + z]]
			msd_pos.append(posiciones[0])
			BLS = simulation.add_topology(topology_type="TopologyC", particle_types=["F1", "exo"],
										  positions=np.array(posiciones))
			graph = BLS.get_graph()
			graph.add_edge(0, 1)

######################## CMC ############################

N = 0

for i in range(0, N, 1):
	x = random.randint(-60, 60)
	y = random.randint(-60, 60)
	z = random.randint(-60, 60)

	# posición del graph de cada hebra del bloque de celulosa
	posiciones = np.array([
		[x + 0., 0. + y, 0. + z],
		[x + 1., 0. + y, 0. + z],
		[x + 2., 0. + y, 0. + z],
		[x + 3., 0. + y, 0. + z],
		[x + 4., 0. + y, 0. + z],
		[x + 5., 0. + y, 0. + z],
		[x + 6., 0. + y, 0. + z],
		[x + 7., 0. + y, 0. + z],
		[x + 8., 0. + y, 0. + z],
		[x + 9., 0. + y, 0. + z],
		[x + 10., 0. + y, 0. + z],
		[x + 11., 0. + y, 0. + z],
		[x + 12., 0. + y, 0. + z],
		[x + 13., 0. + y, 0. + z],
		[x + 14., 0. + y, 0. + z],
		[x + 15., 0. + y, 0. + z],
		[x + 16., 0. + y, 0. + z],
		[x + 17., 0. + y, 0. + z],
		[x + 18., 0. + y, 0. + z],
		[x + 19., 0. + y, 0. + z],
		[x + 20., 0. + y, 0. + z],
		[x + 21., 0. + y, 0. + z],
		[x + 22., 0. + y, 0. + z],
		[x + 23., 0. + y, 0. + z],
		[x + 24., 0. + y, 0. + z],
		[x + 25., 0. + y, 0. + z],
		[x + 26., 0. + y, 0. + z],
		[x + 27., 0. + y, 0. + z],
		[x + 28., 0. + y, 0. + z],
		[x + 29., 0. + y, 0. + z]
	])

	# conformacion del graph de celulosa
	celulosa = simulation.add_topology("TopologyB",
									   ["headA", "core", "core", "core", "core", "core", "core", "core", "core", "core",
										"core", "core", "core", "core", "core", "core", "core", "core", "core", "core",
										"core", "core", "core", "core", "core", "core", "core", "core", "core",
										"headB"], posiciones)

	# Conectividad del graph
	# graph = BLS.get_graph()
	celulosa.get_graph().add_edge(0, 1)
	celulosa.get_graph().add_edge(1, 2)
	celulosa.get_graph().add_edge(2, 3)
	celulosa.get_graph().add_edge(3, 4)
	celulosa.get_graph().add_edge(4, 5)
	celulosa.get_graph().add_edge(5, 6)
	celulosa.get_graph().add_edge(6, 7)
	celulosa.get_graph().add_edge(7, 8)
	celulosa.get_graph().add_edge(8, 9)
	celulosa.get_graph().add_edge(9, 10)
	celulosa.get_graph().add_edge(10, 11)
	celulosa.get_graph().add_edge(11, 12)
	celulosa.get_graph().add_edge(12, 13)
	celulosa.get_graph().add_edge(13, 14)
	celulosa.get_graph().add_edge(14, 15)
	celulosa.get_graph().add_edge(15, 16)
	celulosa.get_graph().add_edge(16, 17)
	celulosa.get_graph().add_edge(17, 18)
	celulosa.get_graph().add_edge(18, 19)
	celulosa.get_graph().add_edge(19, 20)
	celulosa.get_graph().add_edge(20, 21)
	celulosa.get_graph().add_edge(21, 22)
	celulosa.get_graph().add_edge(22, 23)
	celulosa.get_graph().add_edge(23, 24)
	celulosa.get_graph().add_edge(24, 25)
	celulosa.get_graph().add_edge(25, 26)
	celulosa.get_graph().add_edge(26, 27)
	celulosa.get_graph().add_edge(27, 28)
	celulosa.get_graph().add_edge(28, 29)

################################## Crystalline cellulose ##################################

prohibidas = []
def simular_crystalline_cellulose(cantidad):
	N = cantidad*10
	if N == 10:	#Bloque en el centro
		for i in range(0, N, 2):  # original (15x15) ----- loop para hacer un bloque de celulosa
			x = i
			for i in range(0, 10, 2):
				z = i

				# posición del graph de cada hebra del bloque de celulosa
				posiciones = np.array([
					[-15., 2. + x, 0. + z],
					[-14., 2. + x, 0. + z],
					[-13., 2. + x, 0. + z],
					[-12., 2. + x, 0. + z],
					[-11., 2. + x, 0. + z],
					[-10., 2. + x, 0. + z],
					[-9., 2. + x, 0. + z],
					[-8., 2. + x, 0. + z],
					[-7., 2. + x, 0. + z],
					[-6., 2. + x, 0. + z],
					[-5., 2. + x, 0. + z],
					[-4., 2. + x, 0. + z],
					[-3., 2. + x, 0. + z],
					[-2., 2. + x, 0. + z],
					[-1., 2. + x, 0. + z],
					[0., 2. + x, 0. + z],
					[1., 2. + x, 0. + z],
					[2., 2. + x, 0. + z],
					[3., 2. + x, 0. + z],
					[4., 2. + x, 0. + z],
					[5., 2. + x, 0. + z],
					[6., 2. + x, 0. + z],
					[7., 2. + x, 0. + z],
					[8., 2. + x, 0. + z],
					[9., 2. + x, 0. + z],
					[10., 2. + x, 0. + z],
					[11., 2. + x, 0. + z],
					[12., 2. + x, 0. + z],
					[13., 2. + x, 0. + z],
					[14., 2. + x, 0. + z]
				])
				prohibidas.extend([
					[-15., 2. + x, 0. + z],
					[-14., 2. + x, 0. + z],
					[-13., 2. + x, 0. + z],
					[-12., 2. + x, 0. + z],
					[-11., 2. + x, 0. + z],
					[-10., 2. + x, 0. + z],
					[-9., 2. + x, 0. + z],
					[-8., 2. + x, 0. + z],
					[-7., 2. + x, 0. + z],
					[-6., 2. + x, 0. + z],
					[-5., 2. + x, 0. + z],
					[-4., 2. + x, 0. + z],
					[-3., 2. + x, 0. + z],
					[-2., 2. + x, 0. + z],
					[-1., 2. + x, 0. + z],
					[0., 2. + x, 0. + z],
					[1., 2. + x, 0. + z],
					[2., 2. + x, 0. + z],
					[3., 2. + x, 0. + z],
					[4., 2. + x, 0. + z],
					[5., 2. + x, 0. + z],
					[6., 2. + x, 0. + z],
					[7., 2. + x, 0. + z],
					[8., 2. + x, 0. + z],
					[9., 2. + x, 0. + z],
					[10., 2. + x, 0. + z],
					[11., 2. + x, 0. + z],
					[12., 2. + x, 0. + z],
					[13., 2. + x, 0. + z],
					[14., 2. + x, 0. + z]
				])
				# conformacion del graph de celulosa
				celulosa = simulation.add_topology("TopologyB",
												   ["headA", "core", "core", "core", "core", "core", "core", "core", "core",
													"core", "core", "core", "core", "core", "core", "core", "core", "core",
													"core", "core", "core", "core", "core", "core", "core", "core", "core",
													"core", "core", "headB"], posiciones)

				# Conectividad del graph
				# graph = BLS.get_graph()
				celulosa.get_graph().add_edge(0, 1)
				celulosa.get_graph().add_edge(1, 2)
				celulosa.get_graph().add_edge(2, 3)
				celulosa.get_graph().add_edge(3, 4)
				celulosa.get_graph().add_edge(4, 5)
				celulosa.get_graph().add_edge(5, 6)
				celulosa.get_graph().add_edge(6, 7)
				celulosa.get_graph().add_edge(7, 8)
				celulosa.get_graph().add_edge(8, 9)
				celulosa.get_graph().add_edge(9, 10)
				celulosa.get_graph().add_edge(10, 11)
				celulosa.get_graph().add_edge(11, 12)
				celulosa.get_graph().add_edge(12, 13)
				celulosa.get_graph().add_edge(13, 14)
				celulosa.get_graph().add_edge(14, 15)
				celulosa.get_graph().add_edge(15, 16)
				celulosa.get_graph().add_edge(16, 17)
				celulosa.get_graph().add_edge(17, 18)
				celulosa.get_graph().add_edge(18, 19)
				celulosa.get_graph().add_edge(19, 20)
				celulosa.get_graph().add_edge(20, 21)
				celulosa.get_graph().add_edge(21, 22)
				celulosa.get_graph().add_edge(22, 23)
				celulosa.get_graph().add_edge(23, 24)
				celulosa.get_graph().add_edge(24, 25)
				celulosa.get_graph().add_edge(25, 26)
				celulosa.get_graph().add_edge(26, 27)
				celulosa.get_graph().add_edge(27, 28)
				celulosa.get_graph().add_edge(28, 29)

	elif N == 20: #Offset en X para tener dos bloques separados
		for i in range(0, N, 2):  # original (15x15) ----- loop para hacer un bloque de celulosa
			if i < 10:
				x = i + 20 #Offset en el eje X
			else:
				x = i - 20 #Offset en el eje X
			for i in range(0, 10, 2):
				z = i
				# posición del graph de cada hebra del bloque de celulosa
				posiciones = np.array([
					[-15., 2. + x, 0. + z],
					[-14., 2. + x, 0. + z],
					[-13., 2. + x, 0. + z],
					[-12., 2. + x, 0. + z],
					[-11., 2. + x, 0. + z],
					[-10., 2. + x, 0. + z],
					[-9., 2. + x, 0. + z],
					[-8., 2. + x, 0. + z],
					[-7., 2. + x, 0. + z],
					[-6., 2. + x, 0. + z],
					[-5., 2. + x, 0. + z],
					[-4., 2. + x, 0. + z],
					[-3., 2. + x, 0. + z],
					[-2., 2. + x, 0. + z],
					[-1., 2. + x, 0. + z],
					[0., 2. + x, 0. + z],
					[1., 2. + x, 0. + z],
					[2., 2. + x, 0. + z],
					[3., 2. + x, 0. + z],
					[4., 2. + x, 0. + z],
					[5., 2. + x, 0. + z],
					[6., 2. + x, 0. + z],
					[7., 2. + x, 0. + z],
					[8., 2. + x, 0. + z],
					[9., 2. + x, 0. + z],
					[10., 2. + x, 0. + z],
					[11., 2. + x, 0. + z],
					[12., 2. + x, 0. + z],
					[13., 2. + x, 0. + z],
					[14., 2. + x, 0. + z]
				])
				prohibidas.extend([
					[-15., 2. + x, 0. + z],
					[-14., 2. + x, 0. + z],
					[-13., 2. + x, 0. + z],
					[-12., 2. + x, 0. + z],
					[-11., 2. + x, 0. + z],
					[-10., 2. + x, 0. + z],
					[-9., 2. + x, 0. + z],
					[-8., 2. + x, 0. + z],
					[-7., 2. + x, 0. + z],
					[-6., 2. + x, 0. + z],
					[-5., 2. + x, 0. + z],
					[-4., 2. + x, 0. + z],
					[-3., 2. + x, 0. + z],
					[-2., 2. + x, 0. + z],
					[-1., 2. + x, 0. + z],
					[0., 2. + x, 0. + z],
					[1., 2. + x, 0. + z],
					[2., 2. + x, 0. + z],
					[3., 2. + x, 0. + z],
					[4., 2. + x, 0. + z],
					[5., 2. + x, 0. + z],
					[6., 2. + x, 0. + z],
					[7., 2. + x, 0. + z],
					[8., 2. + x, 0. + z],
					[9., 2. + x, 0. + z],
					[10., 2. + x, 0. + z],
					[11., 2. + x, 0. + z],
					[12., 2. + x, 0. + z],
					[13., 2. + x, 0. + z],
					[14., 2. + x, 0. + z]
				])
				# conformacion del graph de celulosa
				celulosa = simulation.add_topology("TopologyB",
												   ["headA", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "headB"], posiciones)

				# Conectividad del graph
				# graph = BLS.get_graph()
				celulosa.get_graph().add_edge(0, 1)
				celulosa.get_graph().add_edge(1, 2)
				celulosa.get_graph().add_edge(2, 3)
				celulosa.get_graph().add_edge(3, 4)
				celulosa.get_graph().add_edge(4, 5)
				celulosa.get_graph().add_edge(5, 6)
				celulosa.get_graph().add_edge(6, 7)
				celulosa.get_graph().add_edge(7, 8)
				celulosa.get_graph().add_edge(8, 9)
				celulosa.get_graph().add_edge(9, 10)
				celulosa.get_graph().add_edge(10, 11)
				celulosa.get_graph().add_edge(11, 12)
				celulosa.get_graph().add_edge(12, 13)
				celulosa.get_graph().add_edge(13, 14)
				celulosa.get_graph().add_edge(14, 15)
				celulosa.get_graph().add_edge(15, 16)
				celulosa.get_graph().add_edge(16, 17)
				celulosa.get_graph().add_edge(17, 18)
				celulosa.get_graph().add_edge(18, 19)
				celulosa.get_graph().add_edge(19, 20)
				celulosa.get_graph().add_edge(20, 21)
				celulosa.get_graph().add_edge(21, 22)
				celulosa.get_graph().add_edge(22, 23)
				celulosa.get_graph().add_edge(23, 24)
				celulosa.get_graph().add_edge(24, 25)
				celulosa.get_graph().add_edge(25, 26)
				celulosa.get_graph().add_edge(26, 27)
				celulosa.get_graph().add_edge(27, 28)
				celulosa.get_graph().add_edge(28, 29)

	elif N == 40: #Offsets en ejes X e Y
		for i in range(0, N, 2):  # original (15x15) ----- loop para hacer un bloque de celulosa
			if i < 20:
				x = i + 30 #Offset en el eje X
				if i < 10:
					y = 30 #Offset en el eje Y
				else:
					y = - 30 #Offset en el eje Y
			else:
				x = i - 30 #Offset en el eje X
				if i < 30:
					y = 30
				else:
					y = - 30
			for i in range(0, 10, 2):
				z = i
				# posición del graph de cada hebra del bloque de celulosa
				posiciones = np.array([
					[-15. + y, 2. + x, 0. + z],
					[-14. + y, 2. + x, 0. + z],
					[-13. + y, 2. + x, 0. + z],
					[-12. + y, 2. + x, 0. + z],
					[-11. + y, 2. + x, 0. + z],
					[-10. + y, 2. + x, 0. + z],
					[-9. + y, 2. + x, 0. + z],
					[-8. + y, 2. + x, 0. + z],
					[-7. + y, 2. + x, 0. + z],
					[-6. + y, 2. + x, 0. + z],
					[-5. + y, 2. + x, 0. + z],
					[-4. + y, 2. + x, 0. + z],
					[-3. + y, 2. + x, 0. + z],
					[-2. + y, 2. + x, 0. + z],
					[-1. + y, 2. + x, 0. + z],
					[0. + y, 2. + x, 0. + z],
					[1. + y, 2. + x, 0. + z],
					[2. + y, 2. + x, 0. + z],
					[3. + y, 2. + x, 0. + z],
					[4. + y, 2. + x, 0. + z],
					[5. + y, 2. + x, 0. + z],
					[6. + y, 2. + x, 0. + z],
					[7. + y, 2. + x, 0. + z],
					[8. + y, 2. + x, 0. + z],
					[9. + y, 2. + x, 0. + z],
					[10. + y, 2. + x, 0. + z],
					[11. + y, 2. + x, 0. + z],
					[12. + y, 2. + x, 0. + z],
					[13. + y, 2. + x, 0. + z],
					[14. + y, 2. + x, 0. + z]
				])
				prohibidas.extend([
					[-15. + y, 2. + x, 0. + z],
					[-14. + y, 2. + x, 0. + z],
					[-13. + y, 2. + x, 0. + z],
					[-12. + y, 2. + x, 0. + z],
					[-11. + y, 2. + x, 0. + z],
					[-10. + y, 2. + x, 0. + z],
					[-9. + y, 2. + x, 0. + z],
					[-8. + y, 2. + x, 0. + z],
					[-7. + y, 2. + x, 0. + z],
					[-6. + y, 2. + x, 0. + z],
					[-5. + y, 2. + x, 0. + z],
					[-4. + y, 2. + x, 0. + z],
					[-3. + y, 2. + x, 0. + z],
					[-2. + y, 2. + x, 0. + z],
					[-1. + y, 2. + x, 0. + z],
					[0. + y, 2. + x, 0. + z],
					[1. + y, 2. + x, 0. + z],
					[2. + y, 2. + x, 0. + z],
					[3. + y, 2. + x, 0. + z],
					[4. + y, 2. + x, 0. + z],
					[5. + y, 2. + x, 0. + z],
					[6. + y, 2. + x, 0. + z],
					[7. + y, 2. + x, 0. + z],
					[8. + y, 2. + x, 0. + z],
					[9. + y, 2. + x, 0. + z],
					[10. + y, 2. + x, 0. + z],
					[11. + y, 2. + x, 0. + z],
					[12. + y, 2. + x, 0. + z],
					[13. + y, 2. + x, 0. + z],
					[14. + y, 2. + x, 0. + z]
				])

				# conformacion del graph de celulosa
				celulosa = simulation.add_topology("TopologyB",
												   ["headA", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "headB"], posiciones)

				# Conectividad del graph
				# graph = BLS.get_graph()
				celulosa.get_graph().add_edge(0, 1)
				celulosa.get_graph().add_edge(1, 2)
				celulosa.get_graph().add_edge(2, 3)
				celulosa.get_graph().add_edge(3, 4)
				celulosa.get_graph().add_edge(4, 5)
				celulosa.get_graph().add_edge(5, 6)
				celulosa.get_graph().add_edge(6, 7)
				celulosa.get_graph().add_edge(7, 8)
				celulosa.get_graph().add_edge(8, 9)
				celulosa.get_graph().add_edge(9, 10)
				celulosa.get_graph().add_edge(10, 11)
				celulosa.get_graph().add_edge(11, 12)
				celulosa.get_graph().add_edge(12, 13)
				celulosa.get_graph().add_edge(13, 14)
				celulosa.get_graph().add_edge(14, 15)
				celulosa.get_graph().add_edge(15, 16)
				celulosa.get_graph().add_edge(16, 17)
				celulosa.get_graph().add_edge(17, 18)
				celulosa.get_graph().add_edge(18, 19)
				celulosa.get_graph().add_edge(19, 20)
				celulosa.get_graph().add_edge(20, 21)
				celulosa.get_graph().add_edge(21, 22)
				celulosa.get_graph().add_edge(22, 23)
				celulosa.get_graph().add_edge(23, 24)
				celulosa.get_graph().add_edge(24, 25)
				celulosa.get_graph().add_edge(25, 26)
				celulosa.get_graph().add_edge(26, 27)
				celulosa.get_graph().add_edge(27, 28)
				celulosa.get_graph().add_edge(28, 29)

	elif N == 80: #Offset en los 3 ejes
		for i in range(0, N, 2):  # original (15x15) ----- loop para hacer un bloque de celulosa
			if i < 40:
				x = i + 30 #Offset en el eje X
				if i < 20:
					y = 30 #Offset en el eje Y
					if i < 10:
						offset = 30 #Offset en el eje Z
					else:
						offset = -30 #Offset en el eje Z
				else:
					y = - 30 #Offset en el eje Y
					if i < 30:
						offset = 30 #Offset en el eje Z
					else:
						offset = -30 #Offset en el eje Z
			else:
				x = i - 30 #Offset en el eje X
				if i < 60:
					y = 30 #Offset en el eje Y
					if i < 50:
						offset = 30 #Offset en el eje Z
					else:
						offset = -30 #Offset en el eje Z
				else:
					y = - 30 #Offset en el eje Y
					if i < 70:
						offset = 30 #Offset en el eje Z
					else:
						offset = -30 #Offset en el eje Z
			for i in range(0, 10, 2):
				z = i
				# posición del graph de cada hebra del bloque de celulosa
				posiciones = np.array([
					[-15. + y, 2. + x, 0. + z + offset],
					[-14. + y, 2. + x, 0. + z+ offset],
					[-13. + y, 2. + x, 0. + z+ offset],
					[-12. + y, 2. + x, 0. + z+ offset],
					[-11. + y, 2. + x, 0. + z+ offset],
					[-10. + y, 2. + x, 0. + z+ offset],
					[-9. + y, 2. + x, 0. + z+ offset],
					[-8. + y, 2. + x, 0. + z+ offset],
					[-7. + y, 2. + x, 0. + z+ offset],
					[-6. + y, 2. + x, 0. + z+ offset],
					[-5. + y, 2. + x, 0. + z+ offset],
					[-4. + y, 2. + x, 0. + z+ offset],
					[-3. + y, 2. + x, 0. + z+ offset],
					[-2. + y, 2. + x, 0. + z+ offset],
					[-1. + y, 2. + x, 0. + z+ offset],
					[0. + y, 2. + x, 0. + z+ offset],
					[1. + y, 2. + x, 0. + z+ offset],
					[2. + y, 2. + x, 0. + z+ offset],
					[3. + y, 2. + x, 0. + z+ offset],
					[4. + y, 2. + x, 0. + z+ offset],
					[5. + y, 2. + x, 0. + z+ offset],
					[6. + y, 2. + x, 0. + z+ offset],
					[7. + y, 2. + x, 0. + z+ offset],
					[8. + y, 2. + x, 0. + z+ offset],
					[9. + y, 2. + x, 0. + z+ offset],
					[10. + y, 2. + x, 0. + z+ offset],
					[11. + y, 2. + x, 0. + z+ offset],
					[12. + y, 2. + x, 0. + z+ offset],
					[13. + y, 2. + x, 0. + z+ offset],
					[14. + y, 2. + x, 0. + z+ offset]
				])
				prohibidas.extend([[-15. + y, 2. + x, 0. + z + offset],
					[-14. + y, 2. + x, 0. + z+ offset],
					[-13. + y, 2. + x, 0. + z+ offset],
					[-12. + y, 2. + x, 0. + z+ offset],
					[-11. + y, 2. + x, 0. + z+ offset],
					[-10. + y, 2. + x, 0. + z+ offset],
					[-9. + y, 2. + x, 0. + z+ offset],
					[-8. + y, 2. + x, 0. + z+ offset],
					[-7. + y, 2. + x, 0. + z+ offset],
					[-6. + y, 2. + x, 0. + z+ offset],
					[-5. + y, 2. + x, 0. + z+ offset],
					[-4. + y, 2. + x, 0. + z+ offset],
					[-3. + y, 2. + x, 0. + z+ offset],
					[-2. + y, 2. + x, 0. + z+ offset],
					[-1. + y, 2. + x, 0. + z+ offset],
					[0. + y, 2. + x, 0. + z+ offset],
					[1. + y, 2. + x, 0. + z+ offset],
					[2. + y, 2. + x, 0. + z+ offset],
					[3. + y, 2. + x, 0. + z+ offset],
					[4. + y, 2. + x, 0. + z+ offset],
					[5. + y, 2. + x, 0. + z+ offset],
					[6. + y, 2. + x, 0. + z+ offset],
					[7. + y, 2. + x, 0. + z+ offset],
					[8. + y, 2. + x, 0. + z+ offset],
					[9. + y, 2. + x, 0. + z+ offset],
					[10. + y, 2. + x, 0. + z+ offset],
					[11. + y, 2. + x, 0. + z+ offset],
					[12. + y, 2. + x, 0. + z+ offset],
					[13. + y, 2. + x, 0. + z+ offset],
					[14. + y, 2. + x, 0. + z+ offset]])

				# conformacion del graph de celulosa
				celulosa = simulation.add_topology("TopologyB",
												   ["headA", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "core", "core", "core", "core", "core", "core",
													"core",
													"core", "core", "headB"], posiciones)

				# Conectividad del graph
				# graph = BLS.get_graph()
				celulosa.get_graph().add_edge(0, 1)
				celulosa.get_graph().add_edge(1, 2)
				celulosa.get_graph().add_edge(2, 3)
				celulosa.get_graph().add_edge(3, 4)
				celulosa.get_graph().add_edge(4, 5)
				celulosa.get_graph().add_edge(5, 6)
				celulosa.get_graph().add_edge(6, 7)
				celulosa.get_graph().add_edge(7, 8)
				celulosa.get_graph().add_edge(8, 9)
				celulosa.get_graph().add_edge(9, 10)
				celulosa.get_graph().add_edge(10, 11)
				celulosa.get_graph().add_edge(11, 12)
				celulosa.get_graph().add_edge(12, 13)
				celulosa.get_graph().add_edge(13, 14)
				celulosa.get_graph().add_edge(14, 15)
				celulosa.get_graph().add_edge(15, 16)
				celulosa.get_graph().add_edge(16, 17)
				celulosa.get_graph().add_edge(17, 18)
				celulosa.get_graph().add_edge(18, 19)
				celulosa.get_graph().add_edge(19, 20)
				celulosa.get_graph().add_edge(20, 21)
				celulosa.get_graph().add_edge(21, 22)
				celulosa.get_graph().add_edge(22, 23)
				celulosa.get_graph().add_edge(23, 24)
				celulosa.get_graph().add_edge(24, 25)
				celulosa.get_graph().add_edge(25, 26)
				celulosa.get_graph().add_edge(26, 27)
				celulosa.get_graph().add_edge(27, 28)
				celulosa.get_graph().add_edge(28, 29)
	else:
		print("TE CONFUNDISTE CON LA CANTIDAD DE SUSTRATO REY")
	return prohibidas

n_bloques_celulosa = 1 #1/2/4/8 valores validos



### SETUP ###


prohibidas = simular_crystalline_cellulose(n_bloques_celulosa)


def main_simular(tipo = sys.argv[2], enzii = sys.argv[3]):
	if tipo == "BLS":
		simular_BLS(enzii)
	elif tipo == "LIBRES":
		simular_LIBRES(enzii)

N_STEPS = 5000000
simulation.make_checkpoints(stride=N_STEPS, output_directory="checkpoints/", max_n_saves=1)
asd = main_simular()
simulation.record_trajectory(stride=1000)
if os.path.exists(simulation.output_file):
	os.remove(simulation.output_file)
simulation.run(n_steps=N_STEPS, timestep=1e-3)


### POST ###
traj = readdy.Trajectory('output')
traj.convert_to_xyz(particle_radii = diccionario)
times, positions = traj.read_observable_particle_positions()
times = np.array(times) * 1e-3


time, counts = traj.read_observable_number_of_particles()

celob = [item[0] for item in counts]
celobB = [item[1] for item in counts]
celobI = [item[2] for item in counts]
gluc = [item[3] for item in counts]
headb = [item[4] for item in counts]
b = [item[5] for item in counts]
celob_and_gluc = [int(item[0]) + int(item[3]) for item in counts]
end = [item[6] for item in counts]
exo = [item[7] for item in counts]

def data_output():
    tim = list(time)
    cel = list(celob)
    celb = list(celobB)
    celI = list(celobI)
    glu = list(gluc)
    hb = list(headb)
    bb = list(b)
    en = list(end)
    ex = list(exo)
    cortes = []
    for i in range(len(tim)):
        corte = cel[i] + celb[i] + celI[i] + glu[i] + hb[i] + bb[i] - hb[1]
        cortes.append(corte) 
    rows=zip(tim,cel,celb,celI,gluc,hb,bb, cortes, en, ex)
    with open("rows.csv", "w") as f:
            writer = csv.writer(f)
            writer.writerow(["time", "cellobiose", "cellobioseB","cellobioseI","glucose", "HeadB", "B", "cortes", "endo", "exo"])
            for row in rows:
                    writer.writerow(row)

data_output()
