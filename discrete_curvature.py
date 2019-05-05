#!/usr/bin/env python
# coding: utf-8

# # Import

# In[1]:


import numpy as np
import numpy.linalg as npla
import sys


# # Helpers

# In[2]:


def get_angle(x_1, x_2):
    
    norm_1 = npla.norm(x_1, 2)
    norm_2 = npla.norm(x_2, 2)
    
    dot_prod = np.dot(x_1, x_2)
    angle = np.arccos(dot_prod / (norm_1 * norm_2))

#     determ = npla.det(np.hstack((x_1.reshape(-1, 1), x_2.reshape(-1, 1))))
#     if detrm > 0: 
#         angle = pi - angle
    
    return angle


# In[3]:


def is_obtuse(v, t, verts, triangs):
    
    flag = 0
    triang = triangs[t]
    x_0, x_1, x_2 = (0, 0, 0)
    
    if triang[0] == v:
        x_1 = np.array(verts[triang[0]])
        x_2, x_3 = np.array(verts[triang[1]]), np.array(verts[triang[2]])
    
    if triang[1] == v:
        x_1 = np.array(verts[triang[1]])
        x_2, x_3 = np.array(verts[triang[0]]), np.array(verts[triang[2]])
    
    if triang[2] == v:
        x_1 = np.array(verts[triang[2]])
        x_2, x_3 = np.array(verts[triang[0]]), np.array(verts[triang[1]])
    
    theta_p = get_angle(x_1-x_2, x_1-x_3)
    if theta_p > np.pi/2:
        flag = 2
    else:
        theta_q, theta_r = get_angle(x_2-x_1, x_2-x_3), get_angle(x_3-x_1, x_3-x_2)
        if theta_q > np.pi/2 or theta_r > np.pi/2:
            flag = 1
    
    return flag


# In[4]:


def get_area(v, t, verts, triangs):
    
    area = 0
    obtuse_flag = is_obtuse(v, t, verts, triangs)
    
    triang_v = triangs[t].copy()
    
    [v_alpha, v_beta] = np.copy(triang_v[triang_v != v])
    
    x = np.array(verts[v])
    x_alpha = np.array(verts[v_alpha])
    x_beta = np.array(verts[v_beta])
        
    if obtuse_flag == 0:
        
        area = 1/8*((npla.norm(x-x_beta, 2)**2)/np.tan(get_angle(x_alpha-x, x_alpha-x_beta)) + (npla.norm(x-x_alpha, 2)**2)/np.tan(get_angle(x_beta-x, x_beta-x_alpha)))
    
    else:
        area = npla.norm(np.cross(x-x_alpha, x-x_beta), 2)/2
        if obtuse_flag == 2:
            area /= 2
        else:
            area /= 4
    
    return area
    


# In[5]:


def get_A_mix(verts, triangs):
    
    A_mix = [0] * len(verts)

    for v in range(len(verts)):
        vert = v
        for t in t_nghbr[v]:
            A_mix[v] += get_area(v, t, verts, triangs)
    
    return A_mix


# In[6]:


def get_mean_curv_norm(A_mix, verts, triangs):
    
    mean_curv_norm = [0] * len(verts)
        
    for v in range(len(verts)):

        vert = v
        for t in t_nghbr[v]:
            
            triang_v = triangs[t].copy()
            
            [v_alpha, v_beta] = np.copy(triang_v[triang_v != v])
            
            x = np.array(verts[v])
            x_alpha = np.array(verts[v_alpha])
            x_beta = np.array(verts[v_beta])
            
            mean_curv_norm[v] += (x-x_beta)/np.tan(get_angle(x_alpha-x, x_alpha-x_beta)) + (x-x_alpha)/np.tan(get_angle(x_beta-x, x_beta-x_alpha))
        
        mean_curv_norm[v] /= 2*A_mix[v]
    
    return np.array(mean_curv_norm)  


# In[7]:


def get_gauss_curv(A_mix, verts, triangs):
    
    gauss_curv = [0] * len(verts)

    for v in range(len(verts)):

        vert = v

        gauss_curv[v] = 2*np.pi

        for t in t_nghbr[v]:

            triang_v = triangs[t].copy()
            
            [v_alpha, v_beta] = np.copy(triang_v[triang_v != v])

            x = np.array(verts[v])
            x_alpha = np.array(verts[v_alpha])
            x_beta = np.array(verts[v_beta])

            gauss_curv[v] -= get_angle(x-x_alpha, x-x_beta)

        gauss_curv[v] /= A_mix[v]

    return np.array(gauss_curv)  


# In[8]:


def get_princ_curv(mean_curv, gauss_curv):
    
    h_curv = mean_curv
    del_curv = (h_curv**2) - gauss_curv
    del_curv[del_curv<0] = 0
        
    min_princ_curv = h_curv - np.sqrt(del_curv)
    max_princ_curv = h_curv + np.sqrt(del_curv)
    
    return min_princ_curv, max_princ_curv


# In[9]:


def read_off(file):
    file = open(file, 'r')
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    n_verts, n_faces, n_dontknow = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')] for i_vert in range(n_verts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for i_face in range(n_faces)]
    return verts, faces


# # Main

# ## Loader

# In[10]:

def main(file):

	verts, triangs = read_off(file)
	# verts, triangs = np.loadtxt('./surf_meshes/sphere/2wcV.txt', dtype='float'), np.loadtxt('./surf_meshes/sphere/2wcT.txt', dtype='int')-1


	# ## Neighbour Vertices

	# In[11]:


	v_nghbr = []

	for i in range(len(verts)):
	    v_nghbr.append([])

	for i in range(len(triangs)):
	    
	    if not triangs[i][1] in v_nghbr[triangs[i][0]] and triangs[i][1] != triangs[i][0]:
	        v_nghbr[triangs[i][0]].append(triangs[i][1])
	    if not triangs[i][2] in v_nghbr[triangs[i][0]] and triangs[i][2] != triangs[i][0]:
	        v_nghbr[triangs[i][0]].append(triangs[i][2])
	    
	    if not triangs[i][0] in v_nghbr[triangs[i][1]] and triangs[i][0] != triangs[i][1]:
	        v_nghbr[triangs[i][1]].append(triangs[i][0])
	    if not triangs[i][2] in v_nghbr[triangs[i][1]] and triangs[i][2] != triangs[i][1]:
	        v_nghbr[triangs[i][1]].append(triangs[i][2])

	    if not triangs[i][0] in v_nghbr[triangs[i][2]] and triangs[i][0] != triangs[i][2]:
	        v_nghbr[triangs[i][2]].append(triangs[i][0])
	    if not triangs[i][1] in v_nghbr[triangs[i][2]] and triangs[i][1] != triangs[i][2]:
	        v_nghbr[triangs[i][2]].append(triangs[i][1])

	# for i in range(len(v_nghbr)):
	#     v_nghbr[i].sort()


	# ## Neighbour Triangles

	# In[12]:


	t_nghbr = []

	for i in range(len(verts)):
	    t_nghbr.append([])

	for i in range(len(triangs)):
	    
	    if triangs[i][0] == triangs[i][1] or triangs[i][1] == triangs[i][2] or triangs[i][0] == triangs[i][2]:
	        continue
	    
	    t_nghbr[triangs[i][0]].append(i)
	    t_nghbr[triangs[i][1]].append(i)
	    t_nghbr[triangs[i][2]].append(i)

	# for i in range(len(t_nghbr)):
	#     t_nghbr[i].sort()


	# # Calc

	# In[13]:


	A_mix = get_A_mix(verts, triangs)

	mean_curv_norm = get_mean_curv_norm(A_mix, verts, triangs)

	mean_curv = npla.norm(mean_curv_norm, ord=2, axis=1)/2

	gauss_curv = get_gauss_curv(A_mix, verts, triangs)

	min_princ_curv, max_princ_curv = get_princ_curv(mean_curv, gauss_curv)


	# ## Save

	# In[14]:


	np.save('K_H.npy', mean_curv)
	np.save('K_G.npy', gauss_curv)
	np.save('K_1.npy', min_princ_curv)
	np.save('K_2.npy', max_princ_curv)
	np.save('vertices.npy', verts)
	np.save('triangles.npy', triangs)

if __name__ == '__main__':

	file_name = sys.argv[1]	
	main(file_name)