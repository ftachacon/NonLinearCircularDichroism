#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.interpolate import interp1d

pi = math.pi
plt.rc('lines', lw=3, color='k')

# In[2]:


lcp = np.loadtxt('lcpSpectrum.txt')
rcp = np.loadtxt('rcpSpectrum.txt')

orderarr = np.loadtxt('orderlist.txt', dtype=int)


# In[3]:


marr = lcp[:, 0]
lcp = lcp[:, 1:]
rcp = rcp[:, 1:]
print(marr)


# In[4]:



Ndat = len(lcp[:, 0])
Norder = len(lcp[0, :])


# In[5]:


cdarr = np.zeros( (Ndat, Norder) )
for i in range(Ndat):
    for j in range(Norder):
        cdarr[i, j] = (lcp[i ,j] - rcp[i, j]) / (lcp[i, j] + rcp[i, j])


# In[6]:


orderarr


# In[7]:


'''averagearr = np.zeros(Ndat)
for i in range(Ndat):
    for j in [0, 2, 4, 6]:
        averagearr[i] += cdarr[i, j]
    averagearr[i] /= 4
fig = plt.figure(figsize=(20,15))
for i in range(Norder):
    if ((orderarr[i]-1)%4 != 0):
        continue
    plt.plot(marr, cdarr[:, i], 'o--', label='H'+str(orderarr[i]))
plt.plot(marr, averagearr, 'ro-', label='average')
plt.xlabel('M', fontsize=20)
plt.ylabel('Circular Dichroism', fontsize=20)
plt.legend()
plt.show()
plt.savefig('CircularDichroism.png')'''


# In[11]:


egap = [0.13780995815870623,
0.11024796652696495,
0.08268597489522372,
0.0551239832634825,
0.02756199163174128,
5.012169450837847e-17,
0.02756199163174128,
0.0551239832634825,
0.08268597489522375,
0.11024796652696497,
0.13780995815870623,
0.11027316781581559,
0.08274178199362908,
0.055240928669274034,
0.027861167929983475,
0.004500657953900623,
0.02799274782920795,
0.05537380884007886,
0.08287490721142261,
0.1104063791207853,
0.13782365194738902,
0.11029840363463261,
0.08279760514651376,
0.05535770729594939,
0.028157323408775383,
0.006365589638445475,
0.028417131571644303,
0.05562259223994828,
0.0830634632894591,
0.11056460493848573,
0.13809001253596823]
# match same parameters with Alexis
corotindex = [0, 2, 4, 6]
width = 11.
hight = width/1.62

fig     = plt.figure(figsize=(width,hight) )
ax      = fig.add_axes([0.135, 0.133, 0.84, 0.8])#,lw=3)
#plt.rc('lines', lw=3, color='k')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);

averagearr = np.zeros(Ndat)
for i in range(Ndat):
    for j in corotindex:
        averagearr[i] += cdarr[i, j]
    averagearr[i] /= len(corotindex)
sw=8
int0=0
p4,=plt.plot(marr, cdarr[:, corotindex[0]], 'ks-', lw=2,markersize=sw )
p5,=plt.plot(marr, cdarr[:, corotindex[1]], 'ys-', lw=2,markersize=sw )
p6,=plt.plot(marr, cdarr[:, corotindex[2]], 'bs-', lw=4,markersize=sw )
p7,=plt.plot(marr, cdarr[:, corotindex[3]], 'gs-', lw=4,markersize=sw )

p9,=plt.plot(marr, averagearr, 'ro-', lw=5 )


plt.legend([p4, p5, p6, p7,  p9], ['HO9','HO13','HO17','HO21','Average'],loc=(0.21, 0.1),fontsize=21)#,loc='lower left')


plt.xlabel(r'$ \rm M $ ', fontsize=34)
plt.ylabel(r'$ \rm Circular\,\,Dichroism $', fontsize=34)


plt.text(-2.85, 1.1, r'$\rm \,(b) \,\,\, Wilson \,\, Mass \,\, Model\,$', fontsize=31, color="k");

xmin    = marr[0]
xmax    = marr[-1]
ymin    = -1.3
ymax    = 1.3
ymax1    = 1.09
shift0  = 0.1
xdown = -2.
xmiddle = 0.
xup = 2.

plt.ylim( ymin, ymax )
plt.xlim( xmin, xmax )
plt.plot([xup,xup],     [ymin,ymax1],    'r--',  lw=4)
plt.plot([xmiddle,xmiddle], [ymin,ymax1],    'r--',  lw=4)
plt.plot([xdown,xdown], [ymin,ymax1],    'r--',  lw=4)

#plt.text(.2, 1.26, "(d)", fontsize=34, color="k")

plt.text(-2.75, -1.25, r"$C=0$", fontsize=26, color="k")
plt.text(xmiddle-1.5, 0.25, r"$C=+1$", fontsize=26, color="k")
plt.text(xup-1.55, 0.25, r"$C=-1$", fontsize=26, color="k")
plt.text(xup+0.18, -1.25, r"$C=0$", fontsize=26, color="k")
plt.tick_params( labelsize=31 )

ax1     = plt.gca()
ax2     = ax1.twiny()

Mt2New  = np.arange(-max(marr),max(marr),2*max(marr)/11.)
ax2.set_xticks( Mt2New )
Npt2    = len(Mt2New)

f       = interp1d( marr, egap, kind='cubic' )
egap1   = f(Mt2New)
ticks1  = []
print('\nLength egap = ', len(egap1), '   Np  = ', Npt2)
for n in range(0,Npt2):
    if n==0:
        ticks1.append(' ')
    else:
        ticks1.append('%.2f'%(egap1[n]*1.))


print('\negap = ',egap1)

plt.yticks([-1,0,1])
ax2.set_xticklabels(ticks1)
ax2.set_xlim( xmin, xmax )
shift0  = 0.1
hd      = -1.1

plt.fill_between( [xmin, xdown-shift0], [1,1], [hd,hd], color=[0.0,0.9,0.], alpha=0.057 )
plt.fill_between( [xdown-shift0, xdown+shift0], [1,1], [hd,hd], color=[0.8,0.0,0.], alpha=0.057 )
plt.fill_between( [xdown+shift0,xmiddle-shift0], [1,1], [hd,hd], color=[0.0,0.8,0.99], alpha=0.04 )
plt.fill_between( [xmiddle-shift0, xmiddle+shift0], [1,1], [hd,hd], color=[0.8,0.0,0.], alpha=0.057 )
plt.fill_between( [xmiddle+shift0,xup-shift0], [1,1], [hd,hd], color=[0.0,0.8,0.99], alpha=0.04 )
plt.fill_between( [xup-shift0, xup+shift0], [1,1], [hd,hd], color=[0.8,0.0,0.], alpha=0.057 )
plt.fill_between( [xup+shift0, xmax], [1,1], [hd,hd], color=[0.0,0.9,0.], alpha=0.057 )

xmin    =  marr[0]
xmax    =  marr[-1]
plt.tick_params( labelsize=18 )
plt.ylim( ymin, ymax )
plt.xlim( xmin, xmax )

fname='FigureApp1b.pdf'
fileName     = fname
plt.savefig(fileName, dpi = 400)

plt.show()


# In[ ]:




