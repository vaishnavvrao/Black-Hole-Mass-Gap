#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/home/BHMG2/anaconda3/bin/python


# In[2]:


import mesa_reader as mr
import numpy as np 
import matplotlib.pyplot as plt


# In[ ]:

#use only when firstly creating the output file
#file= open('output_g_1.txt', 'w')
#file.write("initial_mass,he_dep_mass,CO_core_mass,BH_mass \n")
#file.close()


# In[3]:


def BH(initial_mass):
#change the path accordingly
    history = mr.MesaData('runs_axion/%s/LOGS/history.data' % (initial_mass)) 
    l = mr.MesaLogDir(log_path='runs_axion/%s/LOGS' % (initial_mass))
    prof=l.profile_data()
    
    star_mass=history.star_mass
    co_mass=history.co_core_mass
    center_he4=history.center_he4
    
    num=len(prof.mass)
    BH_mass=0
    e_total=0 #final BH mass
    for i in range(num-1):
        e_total = e_total+ prof.specific_grav_e[i]*(prof.mass[i]-prof.mass[i+1])*1.989*(10**33)
        #print(e_total)
        if e_total<-10**48 and prof.velocity[i]<prof.vesc[i]: #condition for determining BH mass at Fe core collapse
            BH_mass = prof.mass[i]
            break  
            
    index = np.where(center_he4<10**(-3))
    he_depletion= star_mass[index[0][0]-1]
    co_core_mass= co_mass[index[0][0]-1]
    
    file= open('output_g_1.txt', 'a')
    file.write("%s,%s,%s,%s\n" % (initial_mass, he_depletion, co_core_mass, BH_mass))
    file.close()
    
    return 


# In[ ]:

for i in range(40,90,2):
    BH(i)


