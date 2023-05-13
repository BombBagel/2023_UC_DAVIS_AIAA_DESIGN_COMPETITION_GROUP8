import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from labellines import labelLine, labelLines # "pip install matplotlib-label-lines" to install
import sys
import warnings

'''inputs:
W_min and W_max are weights of plane [lbf] or [slug-ft/s^2]
rho_SL and rho_c are density of air at sea level and cruise [slugs/ft^3]
S_ref is wing area [ft^2]
CL_max, CL_min, and CL_alpha [rad^-1] are values based on cruise condition
V_cruise is cruise speed [ft/s]
c_ is mean aerodynamic chord [ft]
h_c is cruise altitude [ft]'''

def vn_diagram(W_min,W_max,rho_SL,rho_c,S_ref,CL_max,CL_min,CL_alpha,V_cruise,c_,h_c):

    def create_line(lowbound,upperbound,value):
        # takes lower and upper bounds at a certain value and makes a horizontal line (flip outputs if you want a vertical line)
        x = np.linspace(lowbound,upperbound)
        y = np.ones(len(x))
        y[:] = value
        return x,y
        
    def create_linear_curve(x1,x2,y1,y2):
        # creates a linear curve between bounds (x1 and x2) with know values for calculating slope
        m = (y2-y1)/(x2-x1)
        b = y1-m*x1
        x = np.linspace(x1,x2,5000)
        y = m*x+b
        return x,y

    def create_gust_velocities_FAR23(x1,x2,y1,y2,scale):
        # creates the FAR23 Specified Gust Velocities Curves
        x,y = create_linear_curve(x1,x2,y1,y2,5000)
        y_tmp,x_tmp = create_line(0,y2,x2)
        y = np.concatenate((y_tmp,y))
        x = np.concatenate((x_tmp,x))
        return x*scale,y*scale

    def create_gust_line(CL_alpha,U_de,W,S_ref,rho,V_max,opt):
        g = 32.1740 # ft/s^2
        
        V = np.linspace(0,V_max,5000)
        if opt == 1: # FAR 25 (lecture slides, unsure if correctly used)
            n_pos,n_neg = FAR25_gust_eqn(CL_alpha,U_de,W,S_ref,rho,V)
        elif opt ==2: # Raymer Equation
            mu = 2*(W/S_ref)/(rho*c_*g*CL_alpha)
            K = 0.88*mu/(5.3+mu)
            n_pos = 1+rho*K*U_de*V*CL_alpha/(2*W/S_ref)
            n_neg = 1-rho*K*U_de*V*CL_alpha/(2*W/S_ref)
        return V,n_pos,n_neg

    def find_intercept(x1, y1, x2, y2):
        # Define a function that takes an x-value and returns the difference between the y-values of the two curves
        def func(x):
            return np.interp(x, x1, y1) - np.interp(x, x2, y2)

        # Choose two initial x-values that bracket the intercept point
        x_left, x_right = np.min([x1, x2]), np.max([x1, x2])

        # Set a tolerance level for the width of the bracketing range
        tol = 1e-8

        # Use the bisection method to iteratively narrow down the range of x-values that bracket the intercept point
        while x_right - x_left > tol:
            x_mid = (x_left + x_right) / 2
            if func(x_mid) > 0:
                x_right = x_mid
            else:
                x_left = x_mid

        # Return the x-coordinate of the intercept point and the corresponding y-coordinate on either curve
        x_intercept = (x_left + x_right) / 2
        y_intercept = np.interp(x_intercept, x1, y1)  # or np.interp(x_intercept, x2, y2)
        return x_intercept, y_intercept

    def FAR25_gust_eqn(CL_alpha,U_de,W,S_ref,rho,V):
        g = 32.1740 # ft/s^2

        mu = 2*(W/S_ref)/(rho*c_*g*CL_alpha)
        K = 0.88*mu/(5.3+mu)
        n_pos = 1+K*CL_alpha*U_de*V*ftpsectoknot/(498*W/S_ref)
        n_neg = 1-K*CL_alpha*U_de*V*ftpsectoknot/(498*W/S_ref)
        return n_pos,n_neg
        
    # V-N Diagram
    weights = [W_min, W_max]
    g = 32.1740 # ft/s^2
    ftpsectoknot = 1/1.688 # knots/ft-s^-1
    nmin = -1 # FAR requirement
    conditions = ["minimum","maximum"]
    font = 18
    sigma = np.sqrt(rho_c/rho_SL)

    for W,condition in zip(weights,conditions):
        V_pos = np.linspace(0,1000,5000)
        V_neg = V_pos.copy()
        n_pos = rho_SL*(V_pos)**2*CL_max/(2*W/S_ref)
        n_neg = rho_SL*(V_neg)**2*CL_min/(2*W/S_ref)

        Vs1 = np.sqrt(2*(W/S_ref)/(rho_SL*CL_max)) # stall speed at 1g
        Vs_1 = np.sqrt(-2*(W/S_ref)/(rho_SL*CL_min)) # stall speed at -1g
        Vc = V_cruise*sigma # factor to convert from true airspeed (TAS) to equivalent airspeed (EAS)
        Vd = 1.25*Vc # dive speed

        if W > 50000: # check if weight of the airplane is larger than 50000 [lbf]
            nd = 2.5
        else:
            low_bound = 2.1+24000/(W+10000)
            nd = (low_bound+3.8)/2 # average FAR requirements
            
        # create stall curves
        index_pos = n_pos <= nd
        index_neg = n_neg >= nmin
        n_pos_final = n_pos[index_pos]
        V_pos_final  = V_pos[index_pos]
        n_neg_final  = n_neg[index_neg]
        V_neg_final  = V_neg[index_neg]

        # create Va to Vd line
        Va = np.sqrt(nd*2*W/S_ref/(rho_SL*CL_max))
        na = nd
        V_ad,n_ad = create_line(Va,Vd,nd)

        # create Vs_1 to Vc line
        V_s_1c,n_s_1c = create_line(Vs_1,Vc,nmin)
        
        # create Vd line
        n_d,V_d = create_line(0,nd,Vd)

        # create Vc to Vd line
        V_cd,n_cd = create_linear_curve(Vc,Vd,nmin,0)
        nc = nd
        
        # plot curves for Vn diagram without gusts
        plt.figure(figsize=(16, 12))
        plt.plot(V_pos_final, n_pos_final, color='red')
        plt.plot(V_neg_final, n_neg_final, color='red')
        plt.plot(V_ad, n_ad, color='red')
        plt.plot(V_s_1c, n_s_1c, color='red')
        plt.plot(V_d, n_d, color='red')
        plt.plot(V_cd, n_cd, color='red')
        #plt.show()
        #plt.clf()

        
        # gust calculations
        ''' # FAR 23 requirements
        scale = 1000
        upperlimit = 20 # 20000 [ft]
        lowerlimit = 50 # 50000 [ft]
        Ue_dive,h_dive = create_gust_velocities_FAR23(12.5,25,lowerlimit,upperlimit,scale)
        Ue_cruise,h_cruise = create_gust_velocities_FAR23(25,50,lowerlimit,upperlimit,scale)
        Ue_rough,h_rough = create_gust_velocities_FAR23(38,66,lowerlimit,upperlimit,scale)'''
        
        # FAR 25 gust speed requirements
        if (h_c >= 0) or (h_c < 15000):
            m = (44-26)/(15000-50000)
            b = 26-m*50000
        elif (h_c >= 15000) or (h_c < 50000):
            m = (56-44)/(0-15000)
            b = 44-m*15000
        U_de_c = m*h_c+b
        U_de_d = U_de_c*0.5 # FAR25 Gust Dive Speed
        
        # cruise
        nc_pos,nc_neg = FAR25_gust_eqn(CL_alpha,U_de_c,W,S_ref,rho_c,Vc)
        
        # dive
        nd_pos,nd_neg = FAR25_gust_eqn(CL_alpha,U_de_d,W,S_ref,rho_c,Vd)

        # recalculating new bounds for the V-n diagram
        mu = 2*(W/S_ref)/(rho_c*c_*g*CL_alpha)
        K = 0.88*mu/(5.3+mu)
        Vb_lim = Vs1*(1+K*U_de_c*(Vc*ftpsectoknot)*CL_alpha/(498*W/S_ref))
        V_b,n_b_pos,n_b_neg = create_gust_line(CL_alpha,U_de_c,W,S_ref,rho_c,1000,1)
        Vb,nb = find_intercept(V_pos,n_pos,V_b,n_b_pos)
        nb_pos,nb_neg = FAR25_gust_eqn(CL_alpha,U_de_c,W,S_ref,rho_c,Vb)
        
        # create Vb to Vc line
        V_bc_pos,n_bc_pos = create_linear_curve(Vb,Vc,nb_pos,nc_pos)
        V_bc_neg,n_bc_neg = create_linear_curve(Vb,Vc,nb_neg,nc_neg)

        # create Vc to Vd line
        V_cd_pos,n_cd_pos = create_linear_curve(Vc,Vd,nc_pos,nd_pos)
        V_cd_neg,n_cd_neg = create_linear_curve(Vc,Vd,nc_neg,nd_neg)
        
        # create Vd line
        n_d,V_d = create_line(nd_neg,nd_pos,Vd)

        # create other gust lines
        V_gust_1,n_gust_1,n_gust_2 = create_gust_line(CL_alpha,U_de_c,W,S_ref,rho_c,Vb,1) # Vb
        V_gust_2 = V_gust_1.copy()
        V_gust_3,n_gust_3,n_gust_4 = create_gust_line(CL_alpha,U_de_c,W,S_ref,rho_c,Vc,1) # Vc
        V_gust_4 = V_gust_3.copy()
        V_gust_5,n_gust_5,n_gust_6 = create_gust_line(CL_alpha,U_de_d,W,S_ref,rho_c,Vd,1) # Vd
        V_gust_6 = V_gust_5.copy()
        
        '''# creates lines directly between points
        V_gust_1,n_gust_1 = create_linear_curve(0,Vb,1,nb_pos) # Vb 
        V_gust_2,n_gust_2 = create_linear_curve(0,Vb,1,nb_neg)
        V_gust_3,n_gust_3 = create_linear_curve(0,Vc,1,nc_pos) # Vc
        V_gust_4,n_gust_4 = create_linear_curve(0,Vc,1,nc_neg)
        V_gust_5,n_gust_5 = create_linear_curve(0,Vd,1,nd_pos) # Vd
        V_gust_6,n_gust_6 = create_linear_curve(0,Vd,1,nd_neg)'''
        
        # create stall curves
        V_int,n_int = find_intercept(V_pos,n_pos,V_gust_2,n_gust_2)
        index_pos = (n_int <= n_pos) & (n_pos <= nb_pos)
        n_pos_final = n_pos[index_pos]
        V_pos_final = V_pos[index_pos]
        
        # create line from intercept of the stall curve to Vb
        V_gust_bot,n_gust_bot = create_linear_curve(V_int,Vb,n_int,nb_neg)

        # check which load factor is larger (manually)
        V_vals = [Vs1,Vs_1,Va,Vc,Vb,Vd]
        n_vals = [1,-1,na,nc_pos,nb,nd]
        V_names = ["Vs1","Vs-1","Va","Vc","Vb","Vd"]
        figname = "Vn Diagram {} Weight with Maneuver and Gust Envelope".format(condition.title())
        
        # plot curves for Vn diagram with gusts
        plt.plot(V_pos_final, n_pos_final, color='blue')
        plt.plot(V_gust_bot, n_gust_bot, color='blue')
        plt.plot(V_bc_pos, n_bc_pos, color='blue')
        plt.plot(V_bc_neg, n_bc_neg, color='blue')
        plt.plot(V_cd_pos, n_cd_pos, color='blue')
        plt.plot(V_cd_neg, n_cd_neg, color='blue')
        plt.plot(V_d, n_d, color='blue')
        plt.plot(V_gust_1,n_gust_1, color='gray',linestyle='--')
        plt.plot(V_gust_2,n_gust_2, color='gray',linestyle='--')
        plt.plot(V_gust_3,n_gust_3, color='gray',linestyle='--',label="{:.2f} [ft/s]".format(U_de_c))
        plt.plot(V_gust_4,n_gust_4, color='gray',linestyle='--')
        plt.plot(V_gust_5,n_gust_5, color='gray',linestyle='--',label="{:.2f} [ft/s]".format(U_de_d))
        plt.plot(V_gust_6,n_gust_6, color='gray',linestyle='--')
        for i,n in enumerate(V_names):
            if (i%2) == 0:
                plt.text(V_vals[i],n_vals[i]+0.05,n,fontsize=font)
            else:
                plt.text(V_vals[i],n_vals[i]-0.1,n,fontsize=font)
            plt.scatter(V_vals[i],n_vals[i],color="black")  
        plt.title("V-n Diagram for {} Weight at {} [ft]".format(condition.title(),h_c),fontsize=font)
        plt.ylabel("n (load factor)",fontsize=font)
        plt.xlabel("Equivalent Air Speed (EAS) [ft/s]",fontsize=font)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            labelLines(plt.gca().get_lines(), align=True, fontsize=font, ha="right")
            
        # create legend
        handles = [plt.Rectangle((0,0), 1, 1, color='red'), plt.Rectangle((0,0), 1, 1, color='blue'), plt.Rectangle((0,0), 1, 1, color='gray')]
        labels = ['Manuever Envelope','Gust Envelope','Gust Lines']
        plt.legend(handles, labels, loc='upper left', fontsize=font)
        
        plt.grid(True)
        plt.xticks(fontsize=font)
        plt.yticks(fontsize=font)
        plt.xlim(0,400)
        #plt.ylim(-1.5,3.5)
        plt.savefig("Vn-diagram-figures/"+figname+".png",dpi=100,bbox_inches='tight')
        #plt.show()
        plt.clf()
        print('\n')
        print("plotted {} weight V-n diagram".format(condition))
        
        # convert back to EAS and print corresponding load factors
        print("Vs 1g = {:.2f}, Vs -1g = {:.2f}, Va = {:.2f}, Vb = {:.2f}, Vc = {:.2f}, Vd = {:.2f} [ft/s] [EAS]".format(Vs1, Vs_1, Va, Vb, Vc, Vd))
        print("n lg = {:.2f}, ns -1g = {:.2f}, na = {:.2f}, nb = {:.2f}, nc = {:.2f}, nd {:.2f}".format(1.0, -1.0, na, nb, nc_pos, nd))
        
        print('\n')
        # WARNINGS
        '''if Vb < Vb_lim:
            Vc_need = (Vb/Vs1-1)*(498*W/S_ref)/(K*U_de_c*CL_alpha)/ftpsectoknot
            print("WARNING V-n Diagram: Vb is too small (Vb = {:.2f} >\= Vs1[1+K*U_c*V_c*CL_alpha/(498*W/S)] = {:.2f} [ft/s] [EAS]), Vc needs to be <= {:.2f} [EAS], {:.2f} [TAS]".format(Vb,Vb_lim,Vc_need,Vc_need/sigma))
        else:
            print("REQUIREMENT SATISFIED: Vb = {:.2f} >= Vs1[1+K*U_c*V_c*CL_alpha/(498*W/S)] = {:.2f} [EAS]".format(Vb,Vb_lim))'''
        if Vc < Vb:
            print("WARNING V-n Diagram: Vc is too small and/or CL_alpha is to large (Vc = {:.2f} >\= Vb = {:.2f} [ft/s] [EAS], CL_alpha = {:.2f} [rad^-1])".format(Vc,Vb,CL_alpha))
        if Va > Vd:
            print("WARNING V-n Diagram: Va is too large (Va = {:.2f} <\= Vd  = {:.2f} [ft/s] [EAS])".format(Va,Vd))
        if Vc < Vs_1:
            print("WARNING V-n Diagram: Vc is too small and/or CL_min is not negative enough (Vc = {:.2f} >\= Vs-1 = {:.2f} [ft/s] [EAS], CL_min = {:.2f}) ".format(Vc,Vs_1,CL_min))
    

# Example
W_min = 52262 # [lbf] or [slug-ft/s^2]
W_max = 64482
rho_SL = 0.00237 # [slugs/ft^3] at sea level
rho_c = 9.61*10**-4 # [slugs/ft^2] at cruise altitude (28000 ft)
#rho_c = 0.001267 # at 20000 ft
S_ref = 730 # [ft^2]
CL_max = 1.55 # for cruise
CL_min = -0.93 # for cruise
CL_alpha = 5.91 # [rad^-1] for cruise
V_cruise = 464.148 # [ft/s] cruise speed (275 knots)
c_ = 8.396 # [ft] mean chord 
h_c = 28000 # [ft] cruise altitude (28000 ft)
vn_diagram(W_min,W_max,rho_SL,rho_c,S_ref,CL_max,CL_min,CL_alpha,V_cruise,c_,h_c)