import numpy as np

graph1 = np.array([[0,1,0,1],	#YUCCA AT THIRD VERTEX
                   [1,0,0,1],
                   [1,1,0,1],
                   [1,0,0,0]])

graph2 = np.array([[0,1,0,1],	#NOT YUCCA
                   [1,0,0,1],
                   [1,1,0,1],
                   [1,1,1,0]])

graph3 = np.array([[0,1,0,0],	#YUCCA AT FOURTH VERTEX
                   [1,0,1,0],
                   [1,1,1,0],
                   [1,0,1,0]])

graph4 = np.array([[0,1,0,1],	#NOT YUCCA
                   [1,0,1,1],
                   [1,1,0,1],
                   [1,0,0,0]])
def yucca(graph):
  for vertex in range(len(graph)):
   if graph[vertex].sum() == len(graph)-1:
     if(graph.sum(axis=0)[vertex])==0:
       return True
  return False  

