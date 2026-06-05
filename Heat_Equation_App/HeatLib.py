import numpy as np
import scipy.linalg as lng
import numexpr as ne
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def Animation(dt,n,t,x,u):
    if n!=0:
        # --- If time step is very small, skip frames to make animation faster ---
        
        if dt<0.001:   
            t_for_animation = t[::10]  
            u_for_animation = u[::10]
        else:
            t_for_animation=t[::1]
            u_for_animation=u[::1]
        fig, ax = plt.subplots()
        
        def animate(i):
            ax.clear()
            ax.plot(x, u_for_animation[i], color="blue")
            ax.set_ylim([u.min(), u.max()])
            ax.set_xlim([x.min(), x.max()])
            ax.set_xlabel("x", fontsize=12)
            ax.set_title("Temperature", fontsize=14, loc="center")  
            ax.grid(True,linewidth=0.5,color="lightgray")
        ani = FuncAnimation(fig=fig,func=animate,frames=len(t_for_animation),interval=1,repeat=False)
        plt.show()
    else:
        return
    
# --- Central-difference spatial discretization using NumPy slicing. --- 

def CDS(r,n,u): return r*u[n,2:]+(1-2*r)*u[n,1:-1]+r*u[n,:-2]

# --- Definition of Forward Time - CDS function with fixed diffusivity a=1 ---
 
def FTCS(TempFormula,a0,aN,f0,fN,k0,kN,dx,dt,L,T,a=1):
    r=(a*dt)/(dx**2)
    x=np.arange(0,L+dx,dx)
    t=np.arange(0,T+dt,dt)
    u=np.zeros((len(t),len(x)))

    # --- Evaluation of formulas for initial and boundary conditions ---

    u[0,:]=ne.evaluate(TempFormula,local_dict={"x":x,"pi":np.pi})
    f0=ne.evaluate(f0,local_dict={"x":x,"t":t,"pi":np.pi})
    fN=ne.evaluate(fN,local_dict={"x":x,"t":t,"pi":np.pi})

    # --- Program control: Select the boundary-condition type from the coefficients a and k ---

    # --- Heat Coefficients are both zero, Dirichlet boundary conditions: k = 0, so u = f/a---

    if k0==0 and kN==0 and a0!=0 and aN!=0:     
        if r>0.5:
            raise ValueError("Unstable problem: FTCS requires r <= 0.5")
        else:
            u[:,0]=f0/a0
            u[:,-1]=fN/aN
            for n in range(0,len(t)-1):
                u[n+1,1:-1]=CDS(r,n,u)
    
    # --- Temperature Coefficients are zero, Neumann boundary conditions: a = 0, so k*u_i = f ---

    elif a0==0 and aN==0 and k0!=0 and kN!=0:   
        if r>0.5:
            raise ValueError("Unstable problem: FTCS requires r <= 0.5")
        else:
            for n in range(0,len(t)-1):
                u[n+1,1:-1]=CDS(r,n,u)
                u[n+1,0]=2*r*u[n,1]+(1-2*r)*u[n,0]-2*dx*r*(f0/k0)
                u[n+1,-1]=2*r*u[n,-2]+(1-2*r)*u[n,-1]-2*dx*r*(fN/kN)

    # --- Robin boundary conditions: k*u_i = f - a*u ----

    elif k0!=0 and kN!=0:                       
        if r>min([1/(2+(a0/k0)*dx),1/(2+(aN/kN)*dx)]):
            raise ValueError("Unstable problem: FTCS requires r <= 0.5")
        else:
            for n in range(0,len(t)-1):
                u[n+1,1:-1]=CDS(r,n,u)
                u[n+1,0]=2*r*u[n,1]+(1-2*r)*u[n,0]+2*dx*r*(a0/k0)*u[n,0]-2*dx*r*(f0/k0)
                u[n+1,-1]=2*r*u[n,-2]+(1-2*r)*u[n,-1]+2*dx*r*(aN/kN)*u[n,-1]-2*dx*r*(fN/kN)
    else:
        n=0
        raise ValueError("Invalid boundary parameters.")

    plt.close()
    Animation(dt,n,t,x,u)

# --- Crank-Nicolson Scheme function ---

def CrankNicolson(TempFormula,a0,aN,f0,fN,k0,kN,dx,dt,L,T,a=1):
    r=(a*dt)/(dx**2)
    x=np.arange(0,L+dx,dx)
    t=np.arange(0,T+dt,dt)
    u=np.zeros((len(t),len(x)))

    u[0,:]=ne.evaluate(TempFormula,local_dict={"x":x,"pi":np.pi})
    f0=ne.evaluate(f0,local_dict={"x":x,"t":t,"pi":np.pi})
    fN=ne.evaluate(fN,local_dict={"x":x,"t":t,"pi":np.pi})
    
    # Matrices for the interior nodes when Dirichlet boundary conditions are used.

    A1=np.diag((2 + 2*r)*np.ones(len(x)-2))+ np.diag(-r*np.ones(len(x)-3),1) + np.diag(-r*np.ones(len(x)-3),-1)
    B1=np.diag((2 - 2*r)*np.ones(len(x)-2))+ np.diag(r*np.ones(len(x)-3),1) + np.diag(r*np.ones(len(x)-3),-1)
    
    # Boundary values are known and are added separately through b0_1 and b1_1.

    b0_1=np.zeros((len(x)-2))
    b1_1=np.zeros_like(b0_1)

    # Matrices for cases where boundary nodes are included in the linear system. This is needed for Neumann and Robin boundary conditions.

    A=np.diag((2 + 2*r)*np.ones(len(x) ))+ np.diag(-r*np.ones(len(x)-1),1) + np.diag(-r*np.ones(len(x)-1),-1) 
    B=np.diag((2 - 2*r)*np.ones(len(x) ))+ np.diag(r*np.ones(len(x)-1),1) + np.diag(r*np.ones(len(x)-1),-1)

    # For Neumann and Robin conditions, b0 and b1 are full time-dependent boundary arrays. 
    # They are filled once for all time levels and indexed

    b1=np.zeros((len(t),len(x)))
    b0=np.zeros_like(b1)

    if k0==0 and kN==0 and a0!=0 and aN!=0:
        if r<=0:
            raise ValueError("Invalid parameters: r must be positive.")
        else:
            u[:,0]=f0/a0
            u[:,-1]=fN/aN

            # For time-dependent Dirichlet boundary conditions, the boundary contribution must be updated at every time step.

            for n in range(0,len(t)-1):
                b0_1[0]=u[n,0]
                b0_1[-1]=u[n,-1]
                b1_1[0]=u[n+1,0]
                b1_1[-1]=u[n+1,-1]
                C1=B1@u[n,1:-1]+r*b0_1+r*b1_1
                u[n+1,1:-1]=lng.solve(A1,C1)

    elif a0==0 and aN==0 and k0!=0 and kN!=0:
        if r<=0:
            raise ValueError("Invalid parameters: r must be positive.")
        else:
            b1[:,0]=b0[:,0]=-2*r*dx*(f0/k0)
            b1[:,-1]=b0[:,-1]=-2*r*dx*(fN/kN)
            A[0,1]=A[-1,-2]=-2*r
            B[0,1]=B[-1,-2]=2*r
            for n in range(len(t)-1):
                C=B@u[n,:]+b0[n,:]+b1[n+1,:]
                u[n+1,:]=lng.solve(A,C)

    elif k0!=0 and kN!=0:
        if r<=0:
            raise ValueError("Invalid parameters: r must be positive.")
        else:
            b1[:,0]=b0[:,0]=2*r*dx*(f0/k0)
            b1[:,-1]=b0[:,-1]=2*r*dx*(fN/kN)
            A[0,1]=A[-1,-2]=-2*r
            B[0,1]=B[-1,-2]=2*r
            A[0,0]=2+2*r*(1-(a0/k0)*dx)
            A[-1,-1]=2+2*r*(1-(aN/kN)*dx)
            B[0,0]=2-2*r*(1-(a0/k0)*dx)
            B[-1,-1]=2-2*r*(1+(aN/kN)*dx)
            for n in range(len(t)-1):
                C=B@u[n,:]+b0[n,:]+b1[n+1,:]
                u[n+1,:]=lng.solve(A,C)
    else:
        n=0
        raise ValueError("Invalid boundary parameters.")

    plt.close()
    Animation(dt,n,t,x,u)

if __name__ == "__main__":
    Initial_Temp=input("Enter initial temperature:")
    f0=input("Left Boundary: ")
    fN=input("Right Boundary: ")
    FTCS(Initial_Temp,1,1,f0,fN,0,0,0.1,0.001,1,1)
