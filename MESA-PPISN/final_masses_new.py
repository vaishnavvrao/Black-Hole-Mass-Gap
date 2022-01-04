#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/home/BHMG2/anaconda3/bin/python


# In[1]:


import mesa_reader as mr
import numpy as np
import matplotlib.pyplot as plt


# In[ ]:


file= open('final_masses.txt', 'w')
file.write("initial_mass,he_dep_mass,co_core_mass,BH_mass \n")
file.close()


# In[2]:


def BH(initial_mass):
    G=6.67428*10**-8
    msun=1.9892*10**33
    rsun=6.9598*10**10
    history = mr.MesaData('runs/%s/LOGS/history.data' % (initial_mass))
    l = mr.MesaLogDir(log_path='runs/%s/LOGS' % (initial_mass))
    prof=l.profile_data()

    star_mass=history.star_mass
    co_mass=history.co_core_mass
    center_he4=history.center_he4 

    num=len(prof.mass) #prof.mass contains interior mass

    dm=np.zeros(num) # list containing masses of individual layers

     for i in range(num-1):
        if i<num-1:
            dm[i]= prof.mass[i]-prof.mass[i+1]
        else:
            dm[i]=prof.mass[i]

    BH_mass= 0

    m_grav=prof.mass*msun
    r=10**(prof.logR)*rsun
    be=G*dm*m_grav/r

    try: 
        BH_mass = m.prof.mass[be>10**48][0]
    except:
        BH_mass= history.M_below_vesc[len(history.M_below_vesc)-1]


    index = np.where(center_he4<10**(-3))
    he_depletion= star_mass[index[0][0]-1]
    co_core_mass= co_mass[index[0][0]-1]

    file= open('compare.txt', 'a')
    file.write("%s,%s,%s,%s\n" % (initial_mass, he_depletion, co_core_mass, BH_mass))
    file.close()

    return


# In[ ]:


for i in range(40,86,2):
    BH(i)
                                                                                                                  87,9          97%
