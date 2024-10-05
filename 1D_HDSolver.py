'''
Astrophysics Gas Dynamics HW5
Brandon Lazard
bjl2698
'''
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.animation as animation
import numpy as np
import pandas as pd

#Define variable arrays and constants
N = 202 #199 cells with the 2 endpoints being ghost cells
t = 0
t_final = 0.5

#Adiabatic index
gamma = 1.4

#Create new and old arrays for all variables
u_new = [0 for i in range(N)]
u_old = [0 for i in range(N)]
p_new = [0 for i in range(N-1)]
p_old = [0 for i in range(N-1)]
rho_new = [0 for i in range(N-1)]
rho_old = [0 for i in range(N-1)]
rho_bar = [0 for i in range(N-1)]
e_new = [0 for i in range(N-1)]
e_old = [0 for i in range(N-1)]
x_new = [0 for i in range(N)]
x_old = [0 for i in range(N)]
q = [0 for i in range(N-1)]
c_s = [0 for i in range(N-1)]
m_cen = [0 for i in range(N-1)]
m_face = [0 for i in range(N)]

#viscosity coefficients
q_0 = 4
q_1 = 0.5

#Domain [a,b]
a = 0
b = 2

#grid spacing
dx = (b - a) / (N-2)

#Apply Symmetry boundary conditions

#Density
rho_old[0] = 1.0
rho_old[-1] = 0.125
#Position
x_old[0] = a
x_old[-1] = a+200*dx
#Mean Density
rho_bar[0] = 1/2
#Pressure
p_old[0] = 1.0
p_old[-1] = 0.1
#Bulk Velocity
u_old[0] = 0
u_old[-1] = 0
#Sound Speed
c_s[0] = np.sqrt((gamma*p_old[0]) / (rho_old[0]))
c_s[-1] = np.sqrt((gamma*p_old[200]) / (rho_old[200]))
#Artificial Viscosity
q[0] = 0
q[-1] = 0

#Mass
m_cen[0] = rho_old[0]*(a+dx - a)

#Set Initial Conditions Everywhere
for i in range(1, N-1):
    if i <= 75: 
        rho_old[i] = 1.0 #Density
        rho_bar[i] = (0.5)*((1/rho_old[i]))
        p_old[i] = 1.0 #Pressure
        c_s[i] = np.sqrt((gamma*p_old[i]) / (rho_old[i])) #sound speed
        q[i] = (c_s[i] / rho_bar[i])*((q_0*((u_old[i+1] - u_old[i])**2))-(q_1*(u_old[i+1] - u_old[i]))) #Artificial Viscosity
        e_old[i] = (p_old[i] / rho_old[i])/(gamma-1) #Specific Energy
        u_old[i] = 0 #Bulk Velocity
        x_old[i] = a + (i-1)*dx #Position
        x_old[i+1] = a + (i*dx) #Position 2
        m_cen[i] = rho_old[i]*(x_old[i+1]-x_old[i]) #Mass defined at the center
        m_face[i] = (1/2)*(m_cen[i-1]+m_cen[i]) #Mass defined at the face
    else:
        rho_old[i] = 0.125 #Density
        rho_bar[i] = (0.5)*((1/rho_old[i])) #Averge density
        p_old[i] = 0.1 #Pressure
        c_s[i] = np.sqrt((gamma*p_old[i]) / (rho_old[i])) #Sound speed
        q[i] = (c_s[i] / rho_bar[i])*((q_0*((u_old[i+1] - u_old[i])**2))-
        (q_1*(u_old[i+1] - u_old[i]))) #Artificial viscosity
        e_old[i] = (p_old[i] / rho_old[i])/(gamma-1) #Specific Energy
        u_old[i] = 0 #Bulk Velocity
        x_old[i] = a + (i-1)*dx #Position
        x_old[i+1] = a + (i*dx) #Position 2
        m_cen[i] = rho_old[i]*(x_old[i+1]-x_old[i]) #Mass defined at the center
        m_face[i] = (1/2)*(m_cen[i-1]+m_cen[i]) #Mass defined at the face

#Mass boundary condition
m_cen[-1] = m_cen[-2]
m_face[-1] = m_face[-2]
m_face[0] = m_face[1]
#Energy boundary conditions
e_old[0] = e_old[1]
e_old[-1] = e_old[-2]

'''
#Define function to update time step
def update_time(cs, u):
cfl = 0.25
dt = cfl*(dx/(max(cs)+max(u)))
return dt
'''
all_den = []
all_en = []
all_p = []
all_u = []
#dt = update_time(c_s, u_old)
dt = 0.001
t = t + dt

#Update lists
while t <= t_final:

    #Create lists that will contain the variables needed to animate the shock wave
    add_den = []
    add_u = []
    add_en = []
    add_p = []

    #Append the initial old element
    add_den.append(rho_old)
    add_u.append(u_old)
    add_en.append(e_old)
    add_p.append(p_old)


    for i in range(1, N-1):
        #time
        dt_cen = (t+dt) - t
        dt_fac = (1/2)*(2*dt_cen)
        #velocity
        u_new[i] = u_old[i] - (dt_fac/m_face[i])*(p_old[i]-p_old[i-1]+q[i]-q[i-1])
        if i == 200:
            u_new[i+1] = u_new[i]
        else:
            u_new[i+1] = u_old[i+1] - (dt_fac/m_face[i+1])*(p_old[i+1]-p_old[i]+q[i+1]-q[i])
        #print(u_new[i])
        #position
        x_new[i] = x_old[i] + dt_cen*u_new[i]
        x_new[i+1] = x_old[i+1] + dt_cen*u_new[i+1]
        #density (FIX X_NEW[I+1])
        rho_new[i] = (m_cen[i]) / (x_new[i+1] - x_new[i])
        #Mean density
        rho_bar[i] = (1/2)*(1/rho_old[i] + 1/rho_new[i])
        #Sound speed
        c_s[i] = np.sqrt(gamma*p_old[i] / rho_old[i])
        #Specific Energy
        e_new[i] = e_old[i] - (p_old[i]+q[i])*(((1/rho_new[i])) - (1/rho_old[i]))
        #Pressure
        p_new[i] = e_new[i]*rho_new[i]*(gamma-1)
        #print(p_new[i])
        
        #Artifical Viscosity
        if (u_old[i+1] - u_old[i]) / (x_old[i+1]-x_old[i]) <= 0:
            q[i] = (c_s[i] / rho_bar[i])*((q_0*((u_old[i+1] - u_old[i])**2))-(q_1*(u_old[i+1] - u_old[i]))) #Artificial Viscosity
        else:
            q[i] = 0
            #Replace all old variables as the new ones
        u_old[i] = u_new[i]
        x_old[i] = x_new[i]
        rho_old[i] = rho_new[i]
        e_old[i] = e_new[i]
        p_old[i] = p_new[i]
        
        add_den.append(rho_new[i])
        add_u.append(u_new[i])
        add_en.append(e_new[i])
        add_p.append(p_new[i])

        

    all_den.append(add_den)
    all_p.append(add_p)
    all_u.append(add_u)
    all_en.append(add_en)

#Update timestep
#dt = update_time(c_s, u_new)
    t = t + dt #Update time
    print(t)



plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

#plt.plot(x_new[1:200], all_den[1][1:200])
#plt.show()

fig, ax = plt.subplots(2, 2, sharex = True, figsize = (15, 12))

def animate(i):

    #print(i)
    #Clear the plots after each frame
    ax[0][0].clear()
    ax[0][1].clear()
    ax[1][0].clear()
    ax[1][1].clear()

    rho = all_den[i]
    u = all_u[i]
    p = all_p[i]
    e = all_en[i]

    ax[0][0].plot(x_new[1:200], rho[1:200], linewidth = 3)
    ax[0][1].plot(x_new[1:200], p[1:200], linewidth = 3)
    ax[1][0].plot(x_new[1:200], e[1:200], linewidth = 3)
    ax[1][1].plot(x_new[1:200], u[1:200], linewidth = 3)

    ax[0][0].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
    ax[0][1].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
    ax[1][0].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
    ax[1][1].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
    ax[1][0].set_xlabel('$Position$', fontsize = 20)
    ax[1][1].set_xlabel('$Position$', fontsize = 20)
    ax[1][0].set_ylabel('$Energy$', fontsize = 20)
    ax[1][1].set_ylabel('$Velocity$', fontsize = 20)
    ax[0][0].set_ylabel('$Density$', fontsize = 20)
    ax[0][1].set_ylabel('$Pressure$', fontsize = 20)
    ax[0][0].grid()
    ax[0][1].grid()
    ax[1][0].grid()
    ax[1][1].grid()

    ax[0][0].set_ylim(0, 1.2)
    ax[0][1].set_ylim(0, 1.2)
    ax[1][0].set_ylim(1.5, 3)
    ax[1][1].set_ylim(0, 1)

ani = animation.FuncAnimation(fig = fig, func = animate, frames = 350, interval = 100)
#ani.save('test.gif')
plt.show()




'''
#Setup Subplots for Density, Pressure, Energy, Velocity
fig, ax = plt.subplots(2, 2, sharex = True, figsize = (15, 12 ))
ax[0][0].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
ax[0][1].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
ax[1][0].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)
ax[1][1].tick_params(which='major', top = True, direction = 'in', length = 6, width = 2, right = True, labelsize = 13)



ax[0][0].plot(x_new[1:200], rho_new[1:200], linewidth = 3)
ax[0][1].plot(x_new[1:200], p_new[1:200], linewidth = 3)
ax[1][0].plot(x_new[1:200], e_new[1:200], linewidth = 3)
ax[1][1].plot(x_new[1:200], u_new[1:200], linewidth = 3)
ax[1][0].set_xlabel('$Position$', fontsize = 20)
ax[1][1].set_xlabel('$Position$', fontsize = 20)
ax[1][0].set_ylabel('$Energy$', fontsize = 20)
ax[1][1].set_ylabel('$Velocity$', fontsize = 20)
ax[0][0].set_ylabel('$Density$', fontsize = 20)
ax[0][1].set_ylabel('$Pressure$', fontsize = 20)
ax[0][0].grid()
ax[0][1].grid()
ax[1][0].grid()
ax[1][1].grid()

plt.savefig('sodshock.jpg', bbox_inches='tight')
plt.show()

plt.plot(x_new[1:201], rho_new[1:201])
plt.plot(x_new[1:201], p_new[1:201])
plt.plot(x_new[1:201], e_new[1:201])
plt.plot(x_new[1:201], u_new[1:201])
'''
