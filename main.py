import pymfinder as py
import time
import pandas as pd
import random
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
package_root = os.path.join(current_dir, "pymfinder-master/pymfinder/")
if package_root not in sys.path:
    sys.path.insert(0, package_root)

import mfinder.mfinder as cmfinder
from itertools import combinations, permutations
from math import sqrt
from roles import *
from datatypes import *
import numpy as np

def motif_structure(network,
                    motifsize = 3,
                    nrandomizations = 0,
                    usemetropolis = False,
                    allmotifs = False,
                    stoufferIDs = False,
                    weighted = False,
                    fweight = None
                    ):

    if motifsize < 2:
        sys.stderr.write("Error: this is not a valid motif size.\n")
        sys.exit()

    if motifsize > 8:
        sys.stderr.write("Error: this is not a recommended motif size.\n")
        sys.exit()

    if motifsize > 4 and allmotifs:
        sys.stderr.write("Warning: 'allmotifs' will be ignored for this motif size and motif_structure will only register existing motifs in the real network.\n")
        allmotifs=False

    if weighted and nrandomizations > 0:
        sys.stderr.write("Warning: the analysis of weighted motifs won't be performed for the randomized networks, only for the real one. There are different ways to randomize weighted networks, you could define your own and run motif_structure multiple times to find the random distribution of weighted motifs.\n")

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stoufferIDs and motifsize!=3:
        sys.stderr.write("Warning: 'stoufferIDs' can only be true when 'motifsize=3' in unipartite networks.\n")
        stats.stoufferIDs = False
    else:
        stats.stoufferIDs = stoufferIDs

    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("Error: you're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

    if stats.networktype:
        stats.networktype = "unipartite"

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("Error: you're trying to mix two different motif analyses.\n")

    stats.weighted = weighted

    if fweight==None:
        fweight = default_fweight

    # parameterize the analysis
    web.MotifSize = motifsize
    web.NRandomizations = nrandomizations
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    #check if this function has already been run
    if len(stats.motifs) != 0:
        stats.motifs = dict()

    # determine all nodes' role statistics
    if stats.weighted:
        if len(stats.motifs) == 0:
            motif_stats(web,stats, allmotifs)
        web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1
        return weighted_motif_stats(web, stats, fweight, allmotifs)
    else:
        return motif_stats(web, stats, allmotifs)

##############################################################
##############################################################
# General utilities
##############################################################
##############################################################

def read_links(filename):
    inFile = open(filename, 'r')
    links = []
    for i in inFile.readlines():
        l = i.strip().split()
        if len(l) > 3 or len(l) < 2:
            inFile.close()
            sys.stderr.write("Error: there is something peculiar about one of the interactions in your input file.\n")
            sys.exit()
        elif len(l) == 2:
            links += [tuple(l + [1])]
        else:
            links += [tuple(l)]

    inFile.close()

    return links

# turn any type of node label into integers (mfinder is finicky like that)
def relabel_nodes(links,stats, buildon=False):
    node_dict = {}
    for i in range(len(links)):
        try:
            s,t,w = links[i]
            w = float(w)
        except ValueError:
            s,t = links[i]
            w = 1

        try:
            s = int(s)
        except:
            pass

        try:
            t = int(t)
        except:
            pass

        if s not in node_dict:
            node_dict[s] = len(node_dict)+1
        if t not in node_dict:
            node_dict[t] = len(node_dict)+1

        if buildon:
            links[i] = (s, t, stats.links[(s,t)].weight)

        else:
            links[i] = (node_dict[s], node_dict[t], w)
            stats.add_link(link_id = (node_dict[s], node_dict[t]), link_name = (s,t))
            stats.links[(node_dict[s], node_dict[t])].weight = w
            try:
                x = stats.nodes[node_dict[s]]
            except KeyError:
                stats.add_node(node_id = node_dict[s], node_name = s)
            try:
                x = stats.nodes[node_dict[t]]
            except KeyError:
                stats.add_node(node_id = node_dict[t], node_name = t)

    return links
    
def gen_mfinder_network(links):
    edges = cmfinder.intArray(len(links)*3+1)
    for i in range(len(links)):
        try:
            s,t,w = links[i]
            w = int(round(w))
        except ValueError:
            s,t = links[i]
            w = 1

        edges[3*i+1] = s
        edges[3*i+2] = t
        edges[3*i+3] = w

    return edges, len(links)

# populate the network info
def mfinder_network_setup(network):
    if type(network) == type("hello world"):
        # DEBUG: if we want to use a filename we need to run a check here to make sure that the node labels are integers and that there are weights
        # web.Filename = network
        stats = NetworkStats()
        network = read_links(network)
        network = relabel_nodes(network,stats, buildon=False)
        edges, numedges = gen_mfinder_network(network)
        return network, stats, edges, numedges
    elif type(network) == type([1,2,3]):
        stats = NetworkStats()
        network = relabel_nodes(network,stats, buildon=False)
        edges, numedges = gen_mfinder_network(network)
        return network, stats, edges, numedges
    elif type(network) == NetworkStats:
        links = relabel_nodes(network.links.keys(), network, buildon=True)
        edges, numedges = gen_mfinder_network(links)
        return links, network, edges, numedges
    else:
        sys.stderr.write("Error: this is an invalid nework input.\n")
        sys.exit()


def default_fweight(x):
    return sum(x)/len(x)


def confidence_interval(data, confidence=0.75):
    av=np.mean(data)
    m=np.median(data)
    sd=np.std(data)
    n=len(data)
    if n==1:
        return data[0], 0.0, data[0],data[0],data[0]
    n_data=np.sort(data)
    mi=n_data[int(round(n*(1-confidence)))]
    ma=n_data[int(round(n*confidence)-1)]
    return av, sd, m, ma, mi


##############################################################
##############################################################
# Motif generating code
##############################################################
##############################################################

def list_motifs(motifsize):

    motifs = cmfinder.list_motifs(motifsize)

    all_motifs = []
    motif_result = motifs.l
    while (motif_result != None):
        all_motifs.append(motif_result.val)
                
        motif_result = motif_result.next

    return all_motifs

def print_motifs(motifsize,motifID=None,outFile=None,links=False,sep=" "):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if motifID:
        motifs=[x for x in list_motifs(motifsize) if int(x)==motifID]
    else:
        motifs=list_motifs(motifsize)

    if motifs==[]:
        sys.stderr.write("Error: this motif does not exist.\n")
        sys.exit()

    for m in motifs:
        output = sep.join(["%i" % m,
                           ])

        fstream.write(output + '\n')

        if links:
            motif_edges = cmfinder.motif_edges(m,motifsize)
            edge_result = motif_edges.l
            while (edge_result != None):
                edge = cmfinder.get_edge(edge_result.p)
                s = int(edge.s)
                t = int(edge.t)
                output = sep.join(["%i" % s,
                                   "%i" % t,
                                   ])
                fstream.write(output + '\n')
                edge_result = edge_result.next

    if outFile:
        fstream.close()

    return   

##############################################################
##############################################################
# random network code
##############################################################
##############################################################

def random_network(network,
                   usemetropolis = False,
                   ):

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # parameterize the analysis
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    return randomized_network(web)
        
def randomized_network(mfinderi):
    results = cmfinder.random_network(mfinderi)
    
    edges = []
    edge_result = results.l
    while (edge_result != None):
        edge = cmfinder.get_edge(edge_result.p)
        s = int(edge.s)
        t = int(edge.t)
        w = int(edge.weight)
        edges.append((s,t,w))

        edge_result = edge_result.next

    return edges

##############################################################
##############################################################
# Motif structure code
##############################################################
##############################################################

def motif_stats(mfinderi,motif_stats, allmotifs):
    results = cmfinder.motif_structure(mfinderi)

    if results:
        motif_result = results.l
        while (motif_result != None):
            motif = cmfinder.get_motif_result(motif_result.p)

            motif_id = int(motif.id)
            
            if allmotifs or int(motif.real_count)+float(motif.rand_mean)!=0:
                motif_stats.add_motif(motif_id)
                motif_stats.motifs[motif_id].real = int(motif.real_count)
                motif_stats.motifs[motif_id].random_m = float(motif.rand_mean)
                motif_stats.motifs[motif_id].random_sd = float(motif.rand_std_dev)
                motif_stats.motifs[motif_id].real_z = float(motif.real_zscore)
                motif_stats.motifs[motif_id].mean_weight = 0.0
                motif_stats.motifs[motif_id].sd_weight = 0.0
                motif_stats.motifs[motif_id].median_weight = 0.0
                motif_stats.motifs[motif_id].firstq_weight = 0.0
                motif_stats.motifs[motif_id].thirdq_weight = 0.0

            motif_result = motif_result.next

    cmfinder.list64_free_mem(results)

    return motif_stats

##############################################################
##############################################################
# Modified code
##############################################################
##############################################################

def participation_stats(mfinderi, participation, links, allmotifs, fweight):
    results = cmfinder.motif_participation(mfinderi)
    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)

    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)

        am_l = motif.all_members.l
        while (am_l != None):
            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            if id in participation.motifs:
                participation.motifs[id].all_members.add(tuple(py_members))

            for m in py_members:
                try:
                    participation.nodes[m].motifs[id] += 1
                except KeyError:
                    participation.nodes[m].motifs[id] = 1

            am_l = am_l.next
        r_l = r_l.next

    cmfinder.res_tbl_mem_free_single(results)
    return participation

def motif_participation(network,
                        links = False,
                        motifsize = 3,
                        randomize = False,
                        usemetropolis = False,
                        stoufferIDs = False,
                        allmotifs = False,
                        weighted = False,
                        fweight = None
                        ):

    if motifsize < 2:
        sys.stderr.write("Error: this is not a valid motif size.\n")
        sys.exit()

    if motifsize > 8:
        sys.stderr.write("Error: this is not a recommended motif size.\n")
        sys.exit()

    if motifsize > 4 and allmotifs:
        sys.stderr.write("Warning: 'allmotifs' will be ignored for this motif size and motif_participation will only register existing motifs in the real network.\n")
        allmotifs=False

    # do we want to randomize the network first?
    if randomize:
        #This will restart the whole object
        network = random_network(network, usemetropolis = usemetropolis)

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stoufferIDs and motifsize!=3:
        sys.stderr.write("Warning: 'stoufferIDs' can only be true when 'motifsize=3' in unipartite networks.\n")
        stats.stoufferIDs = False
    else:
        stats.stoufferIDs = stoufferIDs

    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("Error: you're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

    if stats.networktype:
        stats.networktype = "unipartite"

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("Warning: you're trying to mix two different motif analyses (weighted and not weighted). Be careful!\n")

    stats.weighted = weighted

    if fweight==None:
        fweight = default_fweight

    # parameterize the analysis
    web.MotifSize = motifsize
    web.Randomize = 0
    web.UseMetropolis = 0
    if len(stats.motifs) == 0:
        web.NRandomizations = 0
        web.UseMetropolis = 0
        motif_stats(web,stats, allmotifs)

    web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1

    #TODO I can also run this inside participation
    if stats.weighted:
        weighted_motif_stats(web,stats,fweight,allmotifs)

    #check if this function has already been run
    if len(stats.nodes[stats.nodes.keys()[0]].motifs) != 0:
        for x in stats.nodes.keys():
            stats.nodes[x].motifs = dict()
        if len(stats.links[stats.links.keys()[0]].motifs) != 0:
            for x in stats.links.keys():
                stats.links[x].motifs = dict()

    return participation_stats(web,stats,links,allmotifs,fweight)


import pandas as pd
import time
from multiprocessing import Pool

def _motif_worker(args):
    # Worker function for parallel randomization
    iteration_idx, input_file, id_to_type_map, motif_kwargs = args
    local_data = []
    
    # Generate a random network with preserved degree distribution
    rand_edges = random_network(input_file, usemetropolis=motif_kwargs.get('usemetropolis', False))
    
    # Analyze motif participation in the randomized network
    rand_results = motif_participation(rand_edges, **motif_kwargs)
    for m_id, m_obj in rand_results.motifs.items():
        if m_obj.real > 0 and hasattr(m_obj, 'all_members'):
            for member_ids in m_obj.all_members:
                # Map internal IDs back to original cell type labels
                node_types = [id_to_type_map.get(nid, "Unknown") for nid in member_ids]
                local_data.append([iteration_idx + 1, m_id, str(member_ids), "-".join(node_types)])
    return local_data


def get_motif(input_path, output_path, **motif_kwargs):
    # Analyze the real observed network and save participation results to CSV
    real_results = motif_participation(input_path, **motif_kwargs)
    data_list = []
    for motif_id in sorted(real_results.motifs.keys()):
        m_obj = real_results.motifs[motif_id]
        if m_obj.real > 0 and hasattr(m_obj, 'all_members'):
            for member_ids in m_obj.all_members:
                # Extract original node names
                node_names = [real_results.nodes[nid].id for nid in member_ids]
                data_list.append([motif_id, str(member_ids), "-".join(node_names)])
    df = pd.DataFrame(data_list, columns=['MotifID', 'NodeIDs', 'Celltypes'])
    df.to_csv(output_path, index=False)
    print("Get motif complete. Results saved to {0}".format(output_path))
    return real_results


def motif_random(input_path, output_path, n_randomizations=100, num_cores=40, **motif_kwargs):
    # Initialize mapping
    real_results = motif_participation(input_path, **motif_kwargs)
    id_to_type_map = {nid: str(node_obj.id) for nid, node_obj in real_results.nodes.items()}
    
    # Prepare tasks
    tasks = [(i, input_path, id_to_type_map, motif_kwargs) for i in range(n_randomizations)]
    
    # Manual Pool Management for compatibility
    pool = Pool(processes=num_cores)
    try:
        nested_results = pool.map(_motif_worker, tasks)
    finally:
        pool.close() 
        pool.join()  
        
    # Flatten and save
    random_data_list = [item for sublist in nested_results for item in sublist]
    df_rand = pd.DataFrame(random_data_list, columns=['Iteration', 'MotifID', 'NodeIDs', 'Celltypes'])
    df_rand.to_csv(output_path, index=False)
    
    print("Motif randomization complete. Results saved to {0}".format(output_path))