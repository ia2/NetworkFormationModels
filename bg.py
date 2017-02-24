# This file simulates the network formation model of Bala and Goyal
# new benefits model
from numpy import *
from random import *
import networkx as nx
import matplotlib.pyplot as plt
from itertools import *



def powerset(iterable):
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
         # check http://stackoverflow.com/questions/1482308/whats-a-good-way-to-combinate-through-a-set        
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
        
def cmpperf1(B, tmpgr, Perf): 
        # evaluate the performance of every strategy of node

        for ix,strategy in enumerate(psl):
            # add the edges associated with the strategy
            #print "+++the tested strategy is", strategy
            nnd=[]
            brow=zeros(n)
            for el in strategy:
            #if r is not 0:    
                #print "el is: ", el       
                if el < node:
                    tmpgr.add_edge(node,el)
                else:
                    tmpgr.add_edge(node,el+1)
            #evaluate    the number of direct neigbors
            mud = len(tmpgr.neighbors(node))

            #print "mud is :", mud
            #evaluate the number of nodes whose information is available. 
            ed=nx.bfs_edges(tmpgr,node)
            mu = 0
            for i in ed:
                mu = mu+1
                nnd.append(i[1])
                #update the benefits for node

            for j in range(0,n):
                if j == node:
                    brow[j] = 1
                else:
                    if j in nnd:
                        brow[j]=alphap*B[node,j]+(1-alphap)
                    else:
                        brow[j]=alpham*B[node,j]                            

            #print "mu is :", mu
            MU[ix]=mu
            MUD[ix]=mud
#            Perf[ix]= brow.sum() -c*mud
            Perf[ix]= brow.sum()/(c*max(mud,1))
            # update last performance
            #lastperf[ix,node]=Perf[ix]
            #print "perf is :", Perf[ix]
            #remove the edges for strategy
            for el in strategy:
            #if r is not 0:    
            #print "el is: ", el       
                if el < node:
                   tmpgr.remove_edge(node,el)
                else:
                    tmpgr.remove_edge(node,el+1)   
        return             



def cmpperf2(B, tmpgr, Perf):
        for ix,strategy in enumerate(psl):
            # add the edges associated with the strategy
            #print "+++the tested strategy is", strategy
            nnd=[]
            brow=zeros(n)
            cost=zeros(n)
            for el in strategy:
            #if r is not 0:    
                #print "el is: ", el       
                if el < node:
                    tmpgr.add_edge(node,el)
                else:
                    tmpgr.add_edge(node,el+1)
                
                if el in nb:
                    #print "yes, its a familiar link"
                    cost[el]=c
                else:
                    #print "No, a new link to be established!!"
                    cost[el]=2*c
                    
                    
            #evaluate    the number of direct neigbors
            mud = len(tmpgr.neighbors(node))

            #print "mud is :", mud
            #evaluate the number of nodes whose information is available. 
            ed=nx.bfs_edges(tmpgr,node)
            mu = 0
            for i in ed:
                mu = mu+1
                nnd.append(i[1])
                #update the benefits for node

            for j in range(0,n):
                if j == node:
                    brow[j] = 1
                else:
                    if j in nnd:
                        brow[j]=alphap*B[node,j]+(1-alphap)
                    else:
                        brow[j]=alpham*B[node,j]                            

            #print "mu is :", mu
            print "cost is: ", cost
            MU[ix]=mu
            MUD[ix]=mud
            Perf[ix]= brow.sum() -cost.sum()
#            Perf[ix]= brow.sum()/(c*max(mud,1))
            # update last performance
            #lastperf[ix,node]=Perf[ix]
            #print "perf is :", Perf[ix]
            #remove the edges for strategy
            for el in strategy:
            #if r is not 0:    
            #print "el is: ", el       
                if el < node:
                   tmpgr.remove_edge(node,el)
                else:
                    tmpgr.remove_edge(node,el+1)   
        return             
    

# implement a discounting factor in the graph
def cmpperf3(B, tmpgr, Perf):
        for ix,strategy in enumerate(psl):
            # add the edges associated with the strategy
            #print "+++the tested strategy is", strategy
            nnd=[]
            dben=[] #discounted benefit
            brow=zeros(n)
            cost=zeros(n)
            for el in strategy:
            #if r is not 0:    
                #print "el is: ", el       
                if el < node:
                    tmpgr.add_edge(node,el)
                else:
                    tmpgr.add_edge(node,el+1)
                
                if el in nb:
                    #print "yes, its a familiar link"
                    cost[el]=c
                else:
                    #print "No, a new link to be established!!"
                    cost[el]=2*c
                    
                    
            #evaluate    the number of direct neigbors
            mud = len(tmpgr.neighbors(node))

            #print "mud is :", mud
            #evaluate the number of nodes whose information is available. 
            ed=nx.bfs_edges(tmpgr,node)
            mu = 0
            for i in ed:
                mu = mu+1
                nnd.append(i[1])
                dben.append(len(nx.shortest_path(tmpgr, source=node, target=i[1], weight=None))-1)
                #update the benefits for node
            print " neighbors are", nnd    
            print "distance for discounted benefit", dben
            for j in range(0,n):
                if j == node:
                    brow[j] = 1
                else:
                    if j in nnd:
                        brow[j]=alphap*B[node,j]+(1-alphap)*(delta**dben[nnd.index(j)])
                    else:
                        brow[j]=alpham*B[node,j]                            

            #print "mu is :", mu
            print "cost is: ", cost
            MU[ix]=mu
            MUD[ix]=mud
#            Perf[ix]= brow.sum() - cost.sum()
#            Perf[ix]= brow.sum()/(c*max(mud,1))
            Perf[ix]= brow.sum() - c*mud
            # update last performance
            #lastperf[ix,node]=Perf[ix]
            #print "perf is :", Perf[ix]
            #remove the edges for strategy
            for el in strategy:
            #if r is not 0:    
            #print "el is: ", el       
                if el < node:
                   tmpgr.remove_edge(node,el)
                else:
                    tmpgr.remove_edge(node,el+1)   
        return             
    





def updateG(node, gr, psl, BR):
    tmpnb =     gr[t].neighbors(node)
    # remove the edges associated with the node
    for el in tmpnb:
        gr[t].remove_edge(node,el)
        # add the edges associated with the strategy
    for el in psl[int(BR[node][0])]:
        #if r is not 0:    
        #print "el is: ", el       
        if el < node:
            gr[t].add_edge(node,el)
        else:
            gr[t].add_edge(node,el+1)
    return

def updateB(B,gr,n,node,t,alphap,alpham,tmpgr):
    for i in range(0,n):
        mud = len(tmpgr.neighbors(node))
        ed=nx.bfs_edges(gr[t],i)
        nnd=[]
        mu = 0
        for tmpc in ed:
            mu = mu+1
            nnd.append(tmpc[1])
            #update the benefits for node
            for j in range(0,n):
                if i==j:
                    B[i,j] = 1
                else:
                    if j in nnd:
                        B[i,j]=alphap*B[i,j]+(1-alphap)
                    else:
                        B[i,j]=alpham*B[i,j]  
    return 


n = 4 # number of nodes
c = 0.2
ri = 0.3
alpha = 0.6
alphap=0.8
alpham=0.1
delta = 0.8 # link discount

Timesteps = 1000

gr = []

nset = range(n)

ps = powerset(range(n-1))
psl = list(ps)

# lists containing mu and mud

MU = zeros((len(psl),1))
MUD = zeros((len(psl),1))
Perf = zeros((len(psl),1))
lastperf = zeros((len(psl),n))

B = zeros((n,n)) # benefits matrix
#for i in xrange(1):
for i in xrange(n):
    for j in xrange(n):
        B[i,j]=random()
        #B[i,j]=1
        

lastB = zeros((n,n)) # benefits matrix


BR = zeros((n,1))

simultaneous = 1

for t in range(0,Timesteps):
    
    gr.append(nx.DiGraph())
    # select the initial network randomly
    if t == 0: 
       for node in range(n):
          #print "node is : ", node
           r = randint(0,2**(n-1)-1)
           for el in psl[r]:
               #if r is not 0:    
               #print "el is: ", el       
               if el < node:
                   gr[t].add_edge(node,el)
               else:
                   gr[t].add_edge(node,el+1)
                                   
                            
    else:    
        
        #print "==========This is time step ", t
        #print " The graph edges are ", gr[t-1].edges()
        gr[t]=gr[t-1].copy()
        # select a node randomly to be updated
        if simultaneous:
            for node in xrange(n):
                if random()<=ri:
                    print "node ", node," doesn't update"
                
                else:
                    tmpgr = gr[t-1].copy()
                    nb =     tmpgr.neighbors(node)
                    # remove the edges associated with the node
                    for el in nb:
                        tmpgr.remove_edge(node,el)
                        
                    cmpperf3(B, tmpgr, Perf)
        
                    mx=max(Perf)
                    max_elements=[i for i, j in enumerate(Perf) if j == mx]
                    #print "The maximum elements are ", max_elements
                    #for ff in max_elements:
                    #    print ff
                    br = sample(max_elements,1)
                    print "** For node ", node
                    print "The best response is ", br, "  Strategy is:", psl[br[0]]
                    print "the Brow for the node is: ", B[node,:]
        
                    BR[node] = br
                    # update the graph at time step t with the best response strategy
                    updateG(node, gr, psl, BR)

        
                    print "networks  edges()"
                    print gr[t].edges()
                    #update the B matrix for all nodes
                    updateB(B,gr,n,node, t,alphap,alpham,tmpgr)
                
                
        
        else:
        
            node = randint(0,n-1)
   
            #print "******This is node :", node
        
            tmpgr = gr[t-1].copy()
            nb =     tmpgr.neighbors(node)
            # remove the edges associated with the node
            for el in nb:
                tmpgr.remove_edge(node,el)
            # evaluate the performance of every strategy of node
            
            cmpperf3(B, tmpgr,Perf)
            
            mx=max(Perf)
            max_elements=[i for i, j in enumerate(Perf) if j == mx]
            #print "The maximum elements are ", max_elements
            #for ff in max_elements:
            #    print ff
            br = sample(max_elements,1)
            print "** For node ", node
            print "The best response is ", br, "  Strategy is:", psl[br[0]]
            print "the Brow for the node is: ", B[node,:]
        
            BR[node] = br
            # update the graph at time step t with the best response strategy
            updateG(node, gr, psl, BR)

            #update the B matrix for all nodes
            updateB(B,gr,n,node, t,alphap,alpham,tmpgr)
                        
           
    

# Plot the last 8 network configurations. 
#use interaction mode ion
plt.ion()
fig = plt.figure()
cn = 0
for t in range(-8,0): 
    temp = 420+cn+1
    ax1 = plt.subplot(temp)
    nx.draw_circular(gr[t])
    cn=cn+1
               
    
plt.show()

        
            
            
    
    





