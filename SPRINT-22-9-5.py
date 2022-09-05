#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 10:31:36 2022

@author: liammaher
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:14:02 2022

@author: liammaher
"""

import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from collections import Counter
from itertools import chain
from itertools import product
from itertools import starmap
from functools import partial
chaini = chain.from_iterable
import PySimpleGUI as sg
from PIL import Image



#Function which takes a string as input, splits up the elements by commas and returns a list containing each element
def Convert(string):
    li = list(string.split(","))
    return li



species = 1
ploidy = 2
rows = 1 #number of rows in table

#PySimpleGUI data entry window
layout = [
    [sg.Text('Please enter below a list of taxa names seperated only by commas.')],
    [sg.Text('Taxa:', size =(5, 1)),sg.InputText(size=(50, 1),key=(species,0))],
    [sg.Text('Please enter below a list of the ploidy levels beloging to the above taxa (same order) also seperated only by commas.')],
    [sg.Text('Ploidy:', size =(5, 1)),sg.InputText(size=(50, 1),key=(ploidy,0))], 
    [sg.Text('Alternatively, choose file:', size=(20, 1)), sg.InputText(), sg.FileBrowse(key="-fileIN-")],
    [sg.T('Input Choice:'), sg.Radio('From Text',"input_choice", enable_events=True, default=True, key=("-TEXT-")),
     sg.Radio('From File',"input_choice", enable_events=True, key=("-FILE-"))],
    [sg.T('Node Color:'), sg.Radio('Black',"node_color", enable_events=True, default=True, key=("-IN1-")),
     sg.Radio('Cyan',"node_color", enable_events=True, key=("-IN2-")),
     sg.Radio('Red',"node_color", enable_events=True, key=("-IN3-")),
     sg.Radio('Yellow',"node_color", enable_events=True, key=("-IN4-")),
     sg.Radio('Blue',"node_color", enable_events=True, key=("-IN5-")),
     sg.Radio('Green',"node_color", enable_events=True, key=("-IN6-")),
     ],
    [sg.Text('Custom core network:', size=(20, 1)), sg.InputText(), sg.FileBrowse(key="-filecustom-")],
    [sg.T('Core Network:'), sg.Radio('Optimal',"core_choice", enable_events=True, default=True, key=("-core1-")),
     sg.Radio('Binary',"core_choice", enable_events=True, key=("-core2-")),
     sg.Radio('Prime',"core_choice", enable_events=True, key=("-core3-")),
     sg.Radio('Custom',"core_choice", enable_events=True, key=("-core4-"))
     ],
    [sg.Exit('Finished.', size=(9, 1)), sg.Help('Help')]
    ]
window = sg.Window('Data Entry', layout) 
species_list = str()
ploidy_list = str()
while True:
    event, values = window.read()
    if values["-TEXT-"]==True:
        for x in range(rows):
            species_list = values[(species,x)]
            ploidy_list = values[(ploidy,x)]
        species_list = Convert(species_list)
        ploidy_list = Convert(ploidy_list)
    elif values["-FILE-"]==True:
        file_path = values["-fileIN-"]
        with open(file_path,'r') as f:
            listl=[]
            for line in f:
                strip_lines=line.strip()
                listli = Convert(strip_lines)
                m=listl.append(listli)
            species_list = listl[0]
            ploidy_list = listl[1]
       

#sets the variable n_col to differenct colors depending on the users choice of node color
    if values["-IN1-"]==True:
        n_col = 'black'
    elif values["-IN2-"]==True:
        n_col = 'cyan'
    elif values["-IN3-"]==True:
        n_col = 'red'
    elif values["-IN4-"]==True:
        n_col = 'yellow'
    elif values["-IN5-"]==True:
        n_col = 'blue'
    elif values["-IN6-"]==True:
        n_col = 'green'   
    
# sets the variable core_choice to different values depending on users choice of initialisation
    if values["-core1-"]==True:
        core_choice = 'optimal'
    elif values["-core2-"]==True:
        core_choice = 'binary'
    elif values["-core3-"]==True:
        core_choice = 'prime'
    elif values["-core4-"]==True:
        custom_path = values["-filecustom-"]
        core_choice = 'custom'
        result = []
        with open(custom_path, "r") as fp:
            for i in fp.readlines():
                tmp = i.split(",")
                try:
                    result.append((int(tmp[0]), int(tmp[1]),1))
                except:pass
        customG = nx.MultiDiGraph()
        customG.add_weighted_edges_from(result)
    
 
    
    
    if event == sg.WIN_CLOSED or event == 'Finished.':
        break
    elif event == 'Help':
        sg.Popup("For more information on the theory and results underpinning"+
                 " this program please see The Hybrid Number of a Ploidy Profile"+
                 " written by Dr K. T. Huber and L. J. Maher in 2022."+
                 "\n Contact information: \n k.huber@uea.ac.uk and l.maher@uea.ac.uk.")
    

window.close()  
    

tuples = list(zip(species_list,ploidy_list))
print(tuples)
tuples.sort(key=lambda tup: int(tup[1]), reverse= True)
final_leaves = [x[0] for x in tuples]


# function to take integer as input and return a list of prime factors
def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors


# function which takes an integer as input and returns a connected networkx graph which has
# the input number of beads (sets of parallel arcs) stacked one above the other
def beadLink(num_beads):  
    G = nx.MultiDiGraph()
    i = 0
    while i < num_beads*2:
        if (i % 2) == 0:
            G.add_edge(i,i+1,weight=1.01)
            G.add_edge(i,i+1,weight=5)
            i += 1
        else:
            G.add_edge(i,i+1,weight=1)
            i +=1
    return G

# function which takes an integer input and returns a networkx graph which realises the input using
# the binary representation method (Huber and Maher 2022) by first applying the beadLink function and 
# then operating on this graph.
def strictlySimpleBinary(ploidy_level):
    num = format(ploidy_level, 'b')
    binRep = [int(a) for a in str(num)]
    num_beads = len(binRep)-1
    G = beadLink(num_beads)
    if binRep.count(1) == 1:
        return G

    else:
        num_verts = num_beads*2
        G.remove_edge(0,1)
        G.remove_edge(0,1)
        G.add_edge(0, num_verts+1, weight=1)
        G.add_edge(0, 1, weight=1)
        G.add_edge(num_verts+1, 1, weight=1)
        
        i=1
        for i in range(1,len(binRep)):
            if binRep[-i] == 1:
                low_bead = len(binRep) - i
                hybridvert = 2*low_bead - 1
                break
            else:
                i += 1
                
        G.remove_edge(hybridvert, hybridvert+1)
        G.add_edge(hybridvert, num_verts + 2, weight=1)
        G.add_edge(num_verts + 2, hybridvert + 1, weight=1)
        G.add_edge(num_verts + 1, num_verts + 2, weight=1)
        
        if binRep.count(1) > 2:
            G.remove_edge(num_verts + 1, num_verts + 2)
            subdivis = binRep.count(1) - 2
            G.add_edge(num_verts + 1, num_verts + 3, weight=1)
            G.add_edge(num_verts + 2 + subdivis, num_verts + 2, weight=1)
            i = 1
            for i in range(1,subdivis):
                G.add_edge(num_verts + 2 + i, num_verts + 3 + i, weight=1)
                i += 1

            count = 1
            start_vert = num_verts + 3
            new_hybrid = num_verts + subdivis + 3
            
            binRepLen = len(binRep)
            j = 1
            end_zero_count = 0
            for j in range(1,binRepLen):
                if binRep[-j]==0:
                    end_zero_count += 1
                    j += 1
                else:
                    break
            middleBinRep = binRepLen - 1 - end_zero_count
            for count in range(1,middleBinRep):
                if binRep[count] == 1:
                    hybrid_vert = 2*count - 1
                    if G.has_edge(hybrid_vert,hybrid_vert+1):
                        G.remove_edge(hybrid_vert,hybrid_vert+1)
                        G.add_edge(start_vert,new_hybrid, weight = 1)
                        G.add_edge(hybrid_vert, new_hybrid, weight = 1)
                        G.add_edge(new_hybrid, hybrid_vert+1, weight = 1)
                        count += 1
                        start_vert += 1
                        new_hybrid += 1
                    else:
                        count += 1
                else:
                    count += 1
                    
        return G

# Takes integer input, looks at prime decomposition and 'stitches' together graphs from strictlySimpleBinary
# above one another to realise the input 

def strictlySimplePrime(ploidy_level):
    factors = prime_factors(ploidy_level)
    G = strictlySimpleBinary(factors[0])
    for i in range(1,len(factors)):
        H = strictlySimpleBinary(factors[i])
        H_new = nx.relabel_nodes(H, lambda x:x+len(G))
        F = nx.compose(G,H_new)
        leaves = []
        for node in F.nodes():
            if F.out_degree(node)==0:
                leaves.append(node)
        roots = []
        for node in F.nodes():
            if F.in_degree(node)==0:
                roots.append(node)
        leaf = max(leaves)
        leaf_new = list(F.predecessors(leaf))[0]
        root = min(roots)
        F.remove_edge(leaf_new,leaf)
        F.remove_node(leaf)
        F.add_edge(leaf_new, root, weight = 1)
        mapping = dict(zip(F, range(len(F))))
        G = nx.relabel_nodes(F,mapping)
    return G


# looks at the number of hybrid vertices from both binary and prime factor methods and outputs the 
# graph with fewest
def optimalChoice(ploidy_level):
    A = strictlySimpleBinary(ploidy_level)
    B = strictlySimplePrime(ploidy_level)

    hybrid_list_A = []      
    for node in A.nodes():
        if A.in_degree(node)==2:
            hybrid_list_A.append(node)            
        total_hybrid_number_A = len(hybrid_list_A)
    hybrid_list_B = []      
    for node in B.nodes():
        if B.in_degree(node)==2:
            hybrid_list_B.append(node)            
        
        total_hybrid_number_B = len(hybrid_list_B)
    if total_hybrid_number_A < total_hybrid_number_B:
        G = A
    else:
        G = B    
    return G



# takes a graph which realises a strictly simple integer and adds a tree with leaf size equal to number of
# 1's in the ploidy profile. 
def moveToSimple(ploidy_profile,G):
        ploidy_profile.sort(reverse = True)
        node_list = []
        for node in G.nodes():
            node_list.append(node)
            if G.in_degree(node)==0:
                root = node
                child = list(G.successors(node))[0]
        vertex_addition = max(node_list) +1
        H = nx.MultiDiGraph()
        if len(ploidy_profile) == 2:
            H.add_node(0)
            subdivision_vert = vertex_addition + 1
        else:
            ploidy_one_tree = 2*(len(ploidy_profile) - 1) - 3
            i=0
            subdivision_vert_list = []
            for i in range(ploidy_one_tree):
                if (i % 2) == 0:
                    H.add_edge(i, i+1,weight = 1)
                    H.add_edge(i, i+2, weight = 1)
                    subdivision_vert_list.append(i+2)
                    i+=1
                    subdivision_vert = subdivision_vert_list[-1] + vertex_addition + 1
        H_new = nx.relabel_nodes(H, lambda x:x+vertex_addition)
        F = nx.compose(G,H_new)
        if F[root][child] == {'weight': 1}:
            F.remove_edge(root, child)
            F.remove_edge(root, child)
            F.add_edge(root, subdivision_vert, weight = 1)
            F.add_edge(subdivision_vert, child, weight = 1)
            F.add_edge(subdivision_vert, vertex_addition, weight = 1)
        else:
            F.remove_edge(root, child)
            #F.remove_edge(root, child)
            #F.add_edge(root,child, weight = 1)
            F.add_edge(root, subdivision_vert, weight = 1)
            F.add_edge(subdivision_vert, child, weight = 1)
            F.add_edge(subdivision_vert, vertex_addition, weight = 1)

        return F


def aZero(G,max_h):
    """
    Function which takes a graph made in networkx as input 
    and outputs the resulting graph from the \alpha = 0 case 
    """
    m = max_h
    l = len(G.nodes)
    G.add_edges_from([(m,l),(m,l+1)], weight = 1)
    return G

def aLessThan(G,max_h,alph):
    """
    Function which takes a graph made in networkx as input and 
    outputs the resulting graph from the alpha < m_2 case 
    """
    m = max_h
    a = alph
    l = len(G.nodes)
    G.add_edges_from([(m,l),(m,a),(a,l+1)], weight = 1)
    return G    

def aEqualTo(G,max_h,max2_h):
    """
    Function which takes a graph made in networkx as input and 
    outputs the resulting graph from the alpha = m_2 case 
    """
    m = max_h
    m2 = max2_h
    l = len(G.nodes)
    G.add_edges_from([(m,l),(m,m2),(m2,l+1)], weight = 1)
    return G  

def aGreaterThan(G,max_ploidy_leaf,max2_ploidy_leaf):
    """
    Function which takes a graph made in networkx as input and 
    outputs the resulting graph from the alpha > m_2 case 
    """
    m = max_ploidy_leaf
    m2 = max2_ploidy_leaf
    l = len(G.nodes)
    G.add_edges_from([(m2,m),(m2,l),(m,l+1)], weight = 1)
    return G

simplification_sequence = []

def simp_seq(lst):
    """

    Parameters
    ----------
    lst : List
        Enter the ploidy profile of interest in the form of a list.

    Returns
    -------
    - Each ploidy profile in the simplification sequence 
    paired with the corresponding value of alpha and m_2.
    - The number of hybrid vertices associated 
    with the simplification sequence
    - Lists alpha and m_2 containing their value at 
    each sequence step -- used in the realiseNet function.

    """
    #order lst
    x = sorted(lst, reverse= True)
    simplification_sequence.append(sorted(lst, reverse= True))
    h = 0
    alpha = []
    m_2 = []
    if len(lst) ==1:
        return list(reversed(alpha)), list(reversed(m_2)), list(x)

    else:
    #while the second element of our ploidy profile 
    #is greater than 1 i.e. while ploidy profile \vec{m} is not simple
        while x[1] > 1:
            a = x[0] - x[1]
            # Build lists alpha and m_2
            alpha.append(a)
            m_2.append(x[1])
            #a=0 case doesn't generate a hybrid but we add a 
            #hybrid to the count in both other cases
            if a != 0:
                h += 1
                #replace first element by a
            x[0] = a
            #remove unwanted zeros
            #x = [i for i  in x if i != 0]
            #reorder decreasing
            x = sorted(x, reverse= True)
            c = [i for i  in x if i != 0]
            simplification_sequence.append(c)
        x = [i for i  in x if i != 0]
    return list(reversed(alpha)), list(reversed(m_2)), list(x)


def leaf_choice_equal(G,alpha): 
    """

    Parameters
    ----------
    G : networkx graph
        the input graph that the function uses to determine which vertices are leaves and the number of
        paths from root to those leaves.
    alpha : integer
        the value of alpha at the corresponding step in the simplification sequence for the case alpha = m_2.

    Returns
    -------
    list
        the vertex values of those with the most and second most paths from root to those leaves.

    """
    #assign root and leaf to vertex depending on in/out-degree.
    roots = (v for v, d in G.in_degree() if d == 0)
    leaves = list(v for v, d in G.out_degree() if d == 0)
    #use nx functions to create a list of all paths from root to leaf
    all_paths = partial(nx.all_simple_paths, G)
    leaf_paths = list(chaini(starmap(all_paths, product(roots, leaves))))
    #create new list containing the leaf vertex from each path and print the most common
    new_list = list([sublist[-1] for sublist in leaf_paths])
    max_ploidy = max(set(new_list), key = new_list.count) 
    #Remove all instances of max_ploidy from new_list to obtain the next maximum        
    try:
        while True:
            new_list.remove(max_ploidy)

    except ValueError:
        pass
    max2_ploidy = max(set(new_list), key = new_list.count) 
    return [max_ploidy, max2_ploidy]


def leaf_choice_greater(G,alpha):
    """

    Parameters
    ----------
    G : networkx graph
        the input graph that the function uses to determine which vertices are leaves and the number of
        paths from root to those leaves.
    alpha : integer
        the value of alpha at the corresponding step in the simplification sequence for the case alpha > m_2.

    Returns
    -------
    list
        the vertex values of those with the most and second most paths from root to those leaves.

    """
    #assign root and leaf to vertex depending on in/out-degree.
    roots = (v for v, d in G.in_degree() if d == 0)
    leaves = list(v for v, d in G.out_degree() if d == 0)
    #use nx functions to create a list of all paths from root to leaf
    all_paths = partial(nx.all_simple_paths, G)
    leaf_paths = list(chaini(starmap(all_paths, product(roots, leaves))))
    #create new list containing the leaf vertex from each path and print the most common
    new_list = list([sublist[-1] for sublist in leaf_paths])
    max_ploidy = max(set(new_list), key = new_list.count) 
    #Remove all instances of max_ploidy from new_list to obtain the next maximum        
    try:
        while True:
            new_list.remove(max_ploidy)

    except ValueError:
        pass
    max2_ploidy = max(set(new_list), key = new_list.count) 

    return [max_ploidy, max2_ploidy]


def leaf_choice_zero(G,alpha):
    """

    Parameters
    ----------
    G : networkx graph
        the input graph that the function uses to determine which vertices are leaves and the number of
        paths from root to those leaves.
    alpha : integer
        the value of alpha at the corresponding step in the simplification sequence for the case alpha = 0.

    Returns
    -------
    list
        the vertex values of those with the most and second most paths from root to those leaves.

    """
    #assign root and leaf to vertex depending on in/out-degree.
    roots = (v for v, d in G.in_degree() if d == 0)
    leaves = list(v for v, d in G.out_degree() if d == 0)
    #use nx functions to create a list of all paths from root to leaf
    all_paths = partial(nx.all_simple_paths, G)
    leaf_paths = list(chaini(starmap(all_paths, product(roots, leaves))))
    #create new list containing the leaf vertex from each path and print the most common
    new_list = list([sublist[-1] for sublist in leaf_paths])
    max_ploidy = max(set(new_list), key = new_list.count) 
    #using the same newlist we look for the value which appears alpha times (ploidy level alpha)
    return max_ploidy


def leaf_choice_less(G,alpha):  
    """

    Parameters
    ----------
    G : networkx graph
        the input graph that the function uses to determine which vertices are leaves and the number of
        paths from root to those leaves.
    alpha : integer
        the value of alpha at the corresponding step in the simplification sequence for the case alpha < m_2.

    Returns
    -------
    list
        the vertex values of those with the most and second most paths from root to those leaves.

    """
    #assign root and leaf to vertex depending on in/out-degree.
    roots = (v for v, d in G.in_degree() if d == 0)
    leaves = list(v for v, d in G.out_degree() if d == 0)
    #use nx functions to create a list of all paths from root to leaf
    all_paths = partial(nx.all_simple_paths, G)
    leaf_paths = list(chaini(starmap(all_paths, product(roots, leaves))))
    #create new list containing the leaf vertex from each path and print the most common
    new_list = list([sublist[-1] for sublist in leaf_paths])
    max_ploidy = max(set(new_list), key = new_list.count) 
    #using the same newlist we look for the value which appears alpha times (ploidy level alpha)
    c = Counter(new_list)
    for item in c:
        if c[item] == int(alpha):
            alpha_ploidy = item
            break
    return [max_ploidy, alpha_ploidy]


def SPLINTER(lst):
    """
    

    Parameters
    ----------
    lst : list
        the ploidy profile of interest.

    Returns
    -------
    list
        a networkx graph which realises the input list.

    """
    alpha,m_2,simple_ploidy_profile = simp_seq(lst)        
    ploidy_level = simple_ploidy_profile[0]
    if len(lst) == 1:
        if core_choice == 'optimal':
            G = optimalChoice(ploidy_level)
        elif core_choice == 'binary':
            G = strictlySimpleBinary(ploidy_level)
        elif core_choice == 'prime':
            G = strictlySimplePrime(ploidy_level)
        elif core_choice == 'custom':
            G = customG
    else:
        i=0
        if core_choice == 'optimal':
            G = optimalChoice(ploidy_level)
            if len(simple_ploidy_profile) > 1:
                G = moveToSimple(simple_ploidy_profile,G)
        elif core_choice == 'binary':
            G = strictlySimpleBinary(ploidy_level)
            if len(simple_ploidy_profile) > 1:
                G = moveToSimple(simple_ploidy_profile,G)
        elif core_choice == 'prime':
            G = strictlySimplePrime(ploidy_level)
            if len(simple_ploidy_profile) > 1:
                G = moveToSimple(simple_ploidy_profile,G)
        elif core_choice == 'custom':
            G = customG
      
        
        while(i < len(alpha)):
            #the alpha is equal to 0 case.
            if alpha[i] == 0:
                aZero(G,leaf_choice_zero(G, alpha[i]))
                i += 1
            #the alpha is less than m_2 case.
            elif alpha[i] < m_2[i]:
                output = leaf_choice_less(G, alpha[i])
                aLessThan(G,output[0],output[1])
                i += 1
            #the alpha is equal to m_2 case.
            elif alpha[i] == m_2[i]:
                output = leaf_choice_equal(G, alpha[i])
                aEqualTo(G,output[0],output[1])
                i += 1
            #the alpha is greater than m_2 case.
            elif alpha[i] > m_2[i]:
                output = leaf_choice_greater(G, alpha[i])
                aGreaterThan(G, output[0], output[1])
                i += 1
            else:
                print("error")
                
    hybrid_list = []      
    for node in G.nodes():
        if G.in_degree(node)==2:
            hybrid_list.append(node)            
    total_hybrid_number = len(hybrid_list) 
    return [G,total_hybrid_number]
    
# read data from GUI data entry window as i[0] and i[1] then perform the operations to generate networkx graph
# and species list
G,tot_hybrids = SPLINTER([int(i[1]) for i in tuples])
species = [i[0] for i in tuples]

#Fig drawing material
weights = nx.get_edge_attributes(G,'weight').values()
pos = graphviz_layout(G, prog='dot')
pos2 = nx.spring_layout(G, seed =828282)


# assign the species taxon to correct vertex by equating the number of paths from root to leaf with
# the ploidy levels of the species 
leaves = list(v for v, d in G.out_degree() if d == 0)
roots = (v for v, d in G.in_degree() if d == 0)
leaves = list(v for v, d in G.out_degree() if d == 0)
all_paths = partial(nx.all_simple_paths, G)
leaf_paths = list(chaini(starmap(all_paths, product(roots, leaves))))
new_list = list([sublist[-1] for sublist in leaf_paths])
c = Counter(new_list)
c.most_common()
count_tuples = list(Counter(c).items())
remove_dupes = set(item[1] for item in count_tuples)
distinct_count = [tup for tup in count_tuples if tup[1] in remove_dupes]
from operator import itemgetter
def unique_by_key(elements, key=None):
    if key is None:
        # no key: the whole element must be unique
        key = lambda e: e
    return {key(el): el for el in elements}.values()
distinct_count = list(unique_by_key(count_tuples, key = itemgetter(1)))
L3 = [(x1, y1) for (x1, x2) in tuples for (y1, y2) in count_tuples if int(x2) == int(y2)]
leaf_list = list(dict.fromkeys([i[0] for i in L3]))
vertex_list = list(dict.fromkeys([i[1] for i in L3]))

label= {}

for item in vertex_list:
    pos_leaf = [i for i,x in enumerate(vertex_list) if x == item]
    label[item] = leaf_list[pos_leaf[0]]


# gives arcs contained in a bead curvature so visible in networkx
leftbead = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 1.01]
rightbead = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 5]
rest = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 1]

nx.draw_networkx_nodes(G, pos, node_color = n_col, node_size = 50)
nx.draw_networkx_edges(
    G, pos, edgelist = leftbead,
    connectionstyle="arc3,rad=0.15"  
)
nx.draw_networkx_edges(
    G, pos, edgelist = rightbead,
    connectionstyle="arc3,rad=-0.15"  
)
nx.draw_networkx_edges(
    G, pos, edgelist = rest,
    connectionstyle="arc3,rad=0"
)

pos_higher = {}
y_off = -30  
for k, v in pos.items():
    pos_higher[k] = (v[0], v[1]+y_off)
nx.draw_networkx_labels(G,pos_higher, label)
plt.box(False)
plt.savefig('nx_test.png',bbox_inches='tight')


# resize the figure of the output graph
basewidth = 400
img = Image.open('nx_test.png')
wpercent = (basewidth / float(img.size[0]))
hsize = int((float(img.size[1]) * float(wpercent)))
img = img.resize((basewidth, hsize), Image.ANTIALIAS)
img.save('SPRINT_figure.png')

s = '\n'.join([str(i) for i in simplification_sequence])# pass the string instead of an array.

#pySimpleGUI output window configuration
layout = [

        [sg.Image("SPRINT_figure.png")],
        [sg.Text("The simplification sequence is as follows:")],
        [sg.Text(f"{s}.")],
        [sg.Text(f"The total number of hybrid vertices in this network is {tot_hybrids}.")],
        [sg.B('Save Image'), sg.Exit('Exit', size=(5, 1))]
         ]

window = sg.Window(f"Network that realises {simplification_sequence[0]}.", layout)


while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break      
    else:
        sg.Popup("Image saved in PNG format under file name 'SPRINT_figure.png'.")

window.close()
