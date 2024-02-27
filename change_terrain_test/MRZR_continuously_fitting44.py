from casadi import *
import matplotlib.pyplot as plt
import copy
from numpy import loadtxt
import os.path
import time



def circle(ox, oy, rx,ry):
    theta = np.linspace(0.0, 2*np.pi, 100)
    return ox+rx*np.cos(theta), oy+ry*np.sin(theta)

def CalculateSurrogateFY(FZ, alpha, mu):
    par = np.array([-4.8543,1.4376,1.2205,0.1556,-0.0023,41.3705])
    B= par[0]/mu; 
    C= par[1] 
    D = mu * FZ * par[2]
    E=par[3] 
    Sh = par[4] 
    Sv = par[5]
    X1 = min(max(B * (alpha + Sh), -np.pi/2 + 0.001), np.pi/2 - 0.001)
    F = D * np.sin(C * np.arctan(X1 - E * (X1 - np.arctan(X1)))) + Sv
    return F

def VehicleStatesPredict(states, ctrl):
    la          = 1.5413
    lb          = 1.1740
    M           = 1.3783e+03
    Izz         = 1648
    g           = 9.81
    mu          = 0.6
    KFZF        = 297.9598
    KFZR        = 391.1902
    KFZX        = 160.9084
    KFZYF       = 318.2923
    KFZYR       = 406.5813

    x           = states[0]
    y           = states[1]
    psi         = states[2]
    v           = states[3]
    r           = states[4]
    ux          = states[5]
    sa          = ctrl[0]
    ax          = ctrl[1] 

    zparx = -0.4 * cos(x/15) * sin(y/15)*9/6
    zpary = -0.4 * sin(x/15) * cos(y/15)*9/6
    zparxx = 0.4/15*sin(x/15)*sin(y/15)*9/6
    zparyy = 0.4/15*sin(x/15)*sin(y/15)*9/6
    zparxy = -0.4/15*cos(x/15)*cos(y/15)*9/6


    normal_vec  = np.array([-zparx, -zpary , 1])
    normal_vec  = normal_vec/np.linalg.norm(normal_vec)
    theta       = np.arcsin(normal_vec[0])
    phi         = -np.arcsin(normal_vec[1]/np.arccos(theta))
    RX          = np.array([[1.0,0.0,0.0], [0.0, np.cos(phi), -np.sin(phi)],[0.0, np.sin(phi), np.cos(phi)]])
    RY          = np.array([[np.cos(theta), 0, np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)] ]  )  
    RZ          = np.array([[np.cos(psi),-np.sin(psi), 0],[np.sin(psi),np.cos(psi),0],[0, 0, 1] ]  )   
    R           = RX@ RY@RZ 
    ERX         = np.array([[1,0,np.sin(theta)],[0,np.cos(phi), -np.cos(theta)*np.sin(phi)],[0,np.sin(phi),np.cos(phi)*np.cos(theta)]])

    eb1 = np.array(R[:,0])
    eb2 = np.array(R[:,1])
    eb3 = np.array(R[:,2])

    xps = np.array([1, 0, zparx])
    xpy = np.array([0, 1, zpary])
    xpss = np.array([0, 0, zparxx])
    xpyy = np.array([0, 0, zparyy])
    xpsy = np.array([0, 0, zparxy])
    xpys = np.array([0, 0, zparxy])

    I_matrix = np.array([[np.dot(xps, xps),  np.dot(xps, xpy)],[np.dot(xps, xpy),  np.dot(xpy, xpy)]])
    J_matrix = np.array([[np.dot(xps, eb1),  np.dot(xps, eb2)],[np.dot(eb1, xpy),  np.dot(eb2, xpy)]])
    II_matrix= np.array([[np.dot(xpss, eb3),  np.dot(xpsy, eb3)],[np.dot(xpys, eb3),  np.dot(xpyy, eb3)]])
    b_omega = (np.linalg.pinv(J_matrix)@II_matrix@np.linalg.pinv(I_matrix)@J_matrix)@np.array([ux, v])

    wx          = b_omega[1]
    wy          = -b_omega[0]

    G           = R.T @ np.array([0,0,-g])
    gx          = G[0]
    gy          = G[1]
    gz          = -G[2] + v * wx - ux * wy
    
    FZF         = 2*(KFZF * gz - (ax - r * v - gx) * KFZX)
    FZR         = 2*(KFZR * gz + (ax- r * v - gx) * KFZX)
    FZ          = np.array([FZF, FZR])

    alpha_f     = np.arctan((v + la * r) / (ux + 0.1)) - sa
    alpha_r     = np.arctan((v - lb * r) / (ux + 0.1))

    FY1         = -CalculateSurrogateFY(FZF, -alpha_f, mu)
    FY2         = -CalculateSurrogateFY(FZR, -alpha_r, mu)

    speed       = np.array([ux, v, 0])
    dxdydz      = R @ speed
    dx          = dxdydz[0]
    dy          = dxdydz[1]
    dCardan     = np.linalg.pinv(ERX) @ R @ np.array([wx, wy, r])
    dpsi        = dCardan[2]
    dv          = (FY1 + FY2) / M -  r * ux + gy
    dr          = (FY1 * la - FY2 * lb - 900 * wx * wy) / Izz
    du          = ax
    dstates     = np.array([dx,dy,dpsi,dv,dr,du])

    FZFL_con    = FZF / 2 - ((FY1 + FY2) / M) * KFZYF
    FZFR_con    = FZF / 2 + ((FY1 + FY2) / M) * KFZYF
    FZRL_con    = FZR / 2 - ((FY1 + FY2) / M) * KFZYR
    FZRR_con    = FZR / 2 + ((FY1 + FY2) / M) * KFZYR

    return dstates, np.array([FZFL_con, FZFR_con, FZRL_con, FZRR_con])

def VehicleDynamics3D(states, ctrl, param):
    la          = 1.5413
    lb          = 1.1740
    M           = 1.3783e+03
    Izz         = 1648
    g           = 9.81
    mu          = 0.6
    KFZF        = 297.9598
    KFZR        = 391.1902
    KFZX        = 160.9084
    KFZYF       = 318.2923 * 1.05
    KFZYR       = 406.5813 * 1.05
    rho_ks      = 2.5 
    TS_FZ_circle= 100
    TS_FZ       = 400
    ax_max      = 3.5
    ax_min      = -3.5
    Tire_par    = np.array([-4.8543,1.4376,1.2205,0.1556,-0.0023,41.3705])

    B           = Tire_par[0]/mu; 
    C           = Tire_par[1] 
    E           = Tire_par[3] 
    Sh          = Tire_par[4] 
    Sv          = Tire_par[5]

    x           = states[0]
    y           = states[1]
    psi         = states[2]
    v           = states[3]
    r           = states[4]
    ux          = states[5]
    sa          = states[6]
    sr          = ctrl[0]
    ax          = ctrl[1]


    # zparx = param[0]
    # zpary = param[1]
    # zparxx = param[2]
    # zparyy = param[3]
    # zparxy = param[4]

    # zparx = -0.4 * cos(x/15) * sin(y/15)*9/6
    # zpary = -0.4 * sin(x/15) * cos(y/15)*9/6
    # zparxx = 0.4/15*sin(x/15)*sin(y/15)*9/6
    # zparyy = 0.4/15*sin(x/15)*sin(y/15)*9/6
    # zparxy = -0.4/15*cos(x/15)*cos(y/15)*9/6

    p00         = param[0]
    p10         = param[1]
    p01         = param[2]
    p20         = param[3]
    p11         = param[4]
    p02         = param[5]
    p30         = param[6]
    p21         = param[7]
    p12         = param[8]
    p03         = param[9]
    p40         = param[10]
    p31         = param[11]
    p22         = param[12]
    p13         = param[13]
    p04         = param[14]

    zparx       = (p10 + 2*p20*x + p11*y + 3*p30*(x**2)+2*p21*x*y+p12*(y**2)+4*p40*(x**3)+3*p31*(x**2)*y+2*p22*x*(y**2)+p13*(y**3))
    zpary       = (p01 + 2*p02*y + p11*x + 3*p03*(y**2)+2*p12*x*y+p21*(x**2)+4*p04*(y**3)+3*p13*(y**2)*x+2*p22*y*(x**2)+p31*(x**3))
    zparxx      = 2*p20 + 6*p30*x + 2*p21*y + 12*p40*(x**2) + 6*p31*x*y + 2*p22*(y**2)
    zparyy      = 2*p02 + 6*p03*y + 2*p12*x + 12*p04*(y**2) + 6*p13*y*x + 2*p22*(x**2)
    zparxy      = p11 + 2*p21*x + 2*p12*y + 3*p31*(x**2) + 4*p22*x*y + 3*p13*(y**2)
    zparxy      = p11 + 2*p21*x + 2*p12*y + 3*p31*(x**2) + 4*p22*x*y + 3*p13*(y**2)


    normal_vec  = horzcat(-(zparx), -(zpary) , 1)
    normal_vec  = normal_vec/norm_1(normal_vec)
    theta       = asin(normal_vec[0])
    phi         = -asin(normal_vec[1]/cos(theta))

    RX          = vertcat(horzcat(1.0,0.0,0.0), horzcat(0.0, cos(phi), -sin(phi)),horzcat(0.0, sin(phi), cos(phi)))
    RY          = vertcat( horzcat(cos(theta), 0, sin(theta)), horzcat(0,1,0), horzcat(-sin(theta),0,cos(theta))  )  
    RZ          = vertcat( horzcat(cos(psi),-sin(psi), 0), horzcat(sin(psi),cos(psi),0), horzcat(0, 0, 1) )   
    R           = RX@ RY @RZ 
    ERX         = vertcat( horzcat(1,0,sin(theta)), horzcat(0,cos(phi), -cos(theta)*sin(phi)),horzcat(0,sin(phi),cos(phi)*cos(theta)))
    eb1         = R[:, 0]
    eb2         = R[:, 1]
    eb3         = R[:, 2]
    xps         = vertcat(1, 0, zparx)
    xpy         = vertcat(0, 1, zpary)
    xpss        = vertcat(0, 0, zparxx)
    xpyy        = vertcat(0, 0, zparyy)
    xpsy        = vertcat(0, 0, zparxy)
    xpys        = vertcat(0, 0, zparxy)
    I_matrix    = vertcat(horzcat(dot(xps, xps),  dot(xps, xpy)),horzcat(dot(xps, xpy),  dot(xpy, xpy)))
    J_matrix    = vertcat(horzcat(dot(xps, eb1),  dot(xps, eb2)),horzcat(dot(eb1, xpy),  dot(eb2, xpy)))    
    II_matrix   = vertcat(horzcat(dot(xpss, eb3),  dot(xpsy, eb3)),horzcat(dot(xpys, eb3),  dot(xpyy, eb3)))
    b_omega     = pinv(J_matrix)@II_matrix@pinv(I_matrix)@J_matrix@vertcat(ux, v)

    wx          = b_omega[1]
    wy          = -b_omega[0]

    G           = R.T @ vertcat(0,0,-g)
    gx          = G[0]
    gy          = -g*(cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta))
    gz          = -G[2] + v * wx - ux * wy

    FZF         = 2*(KFZF * gz - (ax - r * v - gx) * KFZX)
    FZR         = 2*(KFZR * gz + (ax- r * v - gx) * KFZX)
    alpha_f     = atan((v + la * r) / (ux + 0.01)) - sa
    alpha_r     = atan((v - lb * r) / (ux + 0.01))    
    X1_f0       = B * (alpha_f + Sh)
    X1_f           = ((( 1 - 1/(1 + exp(-100 * (X1_f0 - 3.1415 / 2)))) + (1/(1 + exp(-100 * (X1_f0 + 3.1415 / 2))))) - 1) * X1_f0 + (1/(1 + exp(-100 * (X1_f0 - 3.1415 / 2)))) * 3.1415 / 2 + (1 - 1/(1 + exp(-100 * (X1_f0 + 3.1415 / 2)))) * -3.1415 / 2
    FY1         = mu * FZF * Tire_par[2] * np.sin(C * atan(X1_f - E * (X1_f - atan(X1_f)))) + Sv
    X1_r0        = B * (alpha_r + Sh)
    X1_r        = ((( 1 - 1/(1 + exp(-100 * (X1_r0 - 3.1415 / 2)))) + (1/(1 + exp(-100 * (X1_r0 + 3.1415 / 2))))) - 1) * X1_r0 + (1/(1 + exp(-100 * (X1_r0 - 3.1415 / 2)))) * 3.1415 / 2 + (1 - 1/(1 + exp(-100 * (X1_r0 + 3.1415 / 2)))) * -3.1415 / 2
    FY2         = mu * FZR * Tire_par[2] * np.sin(C * atan(X1_r - E * (X1_r - atan(X1_r)))) + Sv

    speed       = vertcat(ux, v, 0)
    dxdydz      = R @ speed
    dx          = dxdydz[0]
    dy          = dxdydz[1]
    dCardan     = inv(ERX) @ R @ vertcat(wx, wy, r)
    dpsi        = dCardan[2]
    dv          = (FY1 + FY2) / M -  r * ux + gy
    dr          = (FY1 * la - FY2 * lb) / Izz
    du          = ax
    dsa         = sr
    dstates     = vertcat(dx, dy, dpsi, dv, dr, du, dsa)

    # FZFL_con    = FZF / 2 - ((FY1 + FY2) / M) * KFZYF - TS_FZ
    # FZFR_con    = FZF / 2 + ((FY1 + FY2) / M) * KFZYF - TS_FZ
    # FZRL_con    = FZR / 2 - ((FY1 + FY2) / M) * KFZYR - TS_FZ
    # FZRR_con    = FZR / 2 + ((FY1 + FY2) / M) * KFZYR - TS_FZ

    FZFL_con    = FZF / 2 - (r * ux - gy) * KFZYF - TS_FZ
    FZFR_con    = FZF / 2 + (r * ux - gy) * KFZYF - TS_FZ
    FZRL_con    = FZR / 2 - (r * ux - gy) * KFZYR - TS_FZ
    FZRR_con    = FZR / 2 + (r * ux - gy) * KFZYR - TS_FZ

    FORCE_con1   =  ((ax - r*v - gx) * M ) -  ( ( 2 * KFZF * gz + 2 * KFZR * gz) * mu -  TS_FZ_circle)
    FORCE_con2   = -( ( 2 * KFZF * gz + 2 * KFZR * gz) * mu -  TS_FZ_circle) - ((ax - r*v - gx) * M )
    # FORCE_con   = 0
    ax_max_con  = ax - r*v - gx - ax_max 
    ax_min_con  = ax_min - ax + r * v + gx
    # ks_con         = exp(rho_ks * (- FZFL_con)) + exp(rho_ks * (- FZFR_con)) + exp(rho_ks * (- FZRL_con)) + exp(rho_ks * (- FZRR_con))+ exp(rho_ks * ax_max_con) + exp(rho_ks * ax_min_con)
    ks_con         = exp(rho_ks * ax_max_con) + exp(rho_ks * ax_min_con) #1/rho_ks*( exp(rho_ks * (- FZFL_con)) + exp(rho_ks * (- FZFR_con)) + exp(rho_ks * (- FZRL_con)) + exp(rho_ks * (- FZRR_con)))
    # con = vertcat(ks_con, FORCE_con1, FORCE_con2)
    con = vertcat(ks_con, FORCE_con1, FORCE_con2)
    L = 0.1 * sa**2 + 0.05 * ax**2 + 0.12 * v**2 + 0.3 * sr**2 + 0.02 * (alpha_f**2 + alpha_r**2) + 0.1 * gy**2
    return dstates, con, L


def TransferData(w_opt, OCPParams):
    x_opt = w_opt[1::9]
    y_opt = w_opt[2::9]
    psi_opt = w_opt[3::9]
    v_opt = w_opt[4::9]
    r_opt = w_opt[5::9]
    ux_opt=  w_opt[6::9]
    sa_opt = w_opt[7::9]
    sr_opt = w_opt[8::9]
    ax_opt = w_opt[9::9]
    T =  OCPParams["T"]
    N = OCPParams["Nck"]
    dt = 0.001
    time = np.linspace(0.0, T, N+1)
    detailed_time = np.arange(0, T, dt)
    data_size = detailed_time.shape[0]

    sa_list = numpy.interp(detailed_time, time, sa_opt)
    ax_list = numpy.interp(detailed_time, time, ax_opt)
    sr_list = numpy.interp(detailed_time, time, sr_opt)
    ux_list = numpy.interp(detailed_time, time, ux_opt)
    mat = np.vstack([sa_list, ax_list, ux_list]).T   

    ###################### PLOTTING PREDICTION ##############################
    return mat

def init_states_predict(init_states, controls):
    dt = 1e-3
    data_size = 100
    num_states = init_states.shape[0]
    states_list = np.zeros([data_size, num_states])

    states_list[0,:] = init_states
    cur_states = init_states
    for idx in range(data_size-1):
        ctrl = controls[idx, :2]
        dstates, _ = VehicleStatesPredict(cur_states[:6], ctrl)
        cur_states = cur_states + dt* np.hstack([dstates, 0.0])
        states_list[idx+1,:] = cur_states
        states_list[idx+1,-1] = ctrl[0]
    return states_list[99,:]

def defineSolver(init_states, goal, obslist, OCPParams):
    # Specifty OCP parameters
    nu = 2
    nx = 7
    npar = 15
    OCPParams = {
        "XL":[-2, -2, -pi/2, -2, -pi/2, 1, -pi/9],
        "XU":[102, 102, pi/2, 2, pi/2, 9, pi/9],   
        "CL":[-0.5, -2.5],
        "CU":[0.5, 2.5],
        "XF_S":3.6,
        "OBS_SMX":1.0,
        "OBS_SMY":1.0, 
        "nx":nx,
        "nu":nu,
        "npar":npar,
        "T":3,
        "Nck":10,
        "IntegrationScheme":"Trapezoidal"
    }

    # Declare model variables
    x1 = SX.sym('x1')
    x2 = SX.sym('x2')
    x3 = SX.sym('x3')
    x4 = SX.sym('x4')
    x5 = SX.sym('x5')
    x6 = SX.sym('x6')
    x7 = SX.sym('x7')
    x = vertcat(x1, x2, x3, x4, x5, x6, x7)
    
    u1 = SX.sym('u1')
    u2 = SX.sym('u2')
    u = vertcat(u1, u2)

    par1 = SX.sym('par1')
    par2 = SX.sym('par2')
    par3 = SX.sym('par3')
    par4 = SX.sym('par4')
    par5 = SX.sym('par5')
    par6 = SX.sym('par6')
    par7 = SX.sym('par7')
    par8 = SX.sym('par8')
    par9 = SX.sym('par9')
    par10 = SX.sym('par10')
    par11 = SX.sym('par11')
    par12 = SX.sym('par12')
    par13 = SX.sym('par13')
    par14 = SX.sym('par14')
    par15 = SX.sym('par15')

    parSX = vertcat(par1, par2, par3, par4, par5, par6, par7, par8, par9, par10, par11, par12, par13, par14, par15)


    par = MX.sym("param", npar*(OCPParams["Nck"]+1))

    # param = SX.sym('ter', 5*OCPParams["Nck"])

    xdot, con, L_0 = VehicleDynamics3D(x, u, parSX)

    L = L_0
    f = Function('f', [x, u, parSX], [xdot, con, L], ['x', 'u', 'parSX'], ['xdot', 'con','L'])

    smx = OCPParams["OBS_SMX"]
    smy = OCPParams["OBS_SMY"]

    obslist_sm = copy.copy(obslist)
    obslist_sm[:,2] = obslist_sm[:,2]+smx
    obslist_sm[:,3] = obslist_sm[:,3]+smy

    T = OCPParams["T"]
    N = OCPParams["Nck"]
    nx = OCPParams["nx"]
    nu = OCPParams["nu"]
    XL = OCPParams["XL"]
    XU = OCPParams["XU"]
    CL = OCPParams["CL"]
    CU = OCPParams["CU"]
    dt = T/(N)
    mu = MX.sym('mu'); # slack variable

    w=[]
    w0 = []
    lbw = []
    ubw = []
    J = 0
    g=[]
    lbg = []
    ubg = []

    w.append(mu)
    lbw.append(0)
    ubw.append(inf)
    w0.append(1e-3)


    Xk = MX.sym('X0', nx)
    w.append(Xk)
    if init_states[5] <= XL[5]:
        init_states[5] = XL[5]
    elif init_states[5] >= XU[5]:
        init_states[5] = XU[5]

    if init_states[6] <= XL[6]:
        init_states[6] = XL[6]
    elif init_states[6] >= XU[6]:
        init_states[6] = XU[6]

    mystates = init_states


    lbw.append(mystates)
    ubw.append(mystates)
    w0.append(mystates)
    Uk = MX.sym('U_0', nu)
    w.append(Uk)
    lbw.append(CL)
    ubw.append(CU)
    w0.append([0, 0])



    [fj, conj, qj] = f(Xk, Uk, par[:npar])
    # [fj, conj, qj] = f(Xk, Uk)

    J = J + qj*dt
    for k in range(N):
        Xkpl = MX.sym('X_' + str(k+1), nx)
        w.append(Xkpl)
        lbw.append(XL)
        ubw.append(XU)
        w0.append(mystates)
        Ukpl = MX.sym('U_' + str(k+1), nu)
        w.append(Ukpl)
        lbw.append(CL)
        ubw.append(CU)
        w0.append([1e-3, 0])
        # [fjpl, conjpl, qjpl] = f(Xkpl, Ukpl)
        [fjpl, conjpl, qjpl] = f(Xkpl, Ukpl, par[(k+1)*npar:(k+2)*npar])
        J = J + qjpl*dt

        if (OCPParams["IntegrationScheme"] == "Trapezoidal"):
            g.append(Xkpl - (fjpl+fj)/2 *dt -Xk)
        elif (OCPParams["IntegrationScheme"] == "BkwEuler"):
            g.append(Xkpl - (fjpl)/2 *dt -Xk)
        else:
            raise Exception("No integration scheme available")

        lbg.append(np.zeros(nx))
        ubg.append(np.zeros(nx))

        g.append(conj)
        lbg.append([0.0, -inf, -inf])
        ubg.append([1.0, 0.0, 0.0])

        for obs_index in range(obslist_sm.shape[0]):
            g.append((   (Xkpl[0]-obslist_sm[obs_index,0])/obslist_sm[obs_index,2]   )**2+( (Xkpl[1]-obslist_sm[obs_index,1])/obslist_sm[obs_index,3] )**2)
            lbg.append(1.0)
            ubg.append(inf)

        Xk = Xkpl
        Uk = Ukpl
        fj = fjpl
        qj = qjpl
        conj = conjpl


    J = J+(1e3)*mu**2 + 125.0*(  ((Xk[0]-goal[0])**2+(Xk[1] - goal[1])**2)/((mystates[0]-goal[0])**2+(mystates[1] - goal[1])**2)    )
    w = vertcat(*w)
    g = vertcat(*g)
    prob = {'f': J, 'x': w, 'g': g, 'p':par}
    # opts = {'ipopt.print_level': 0,
    #         'ipopt.tol' : 4e-2,
    #         'ipopt.warm_start_init_point' : 'yes',
    #         # 'ipopt.warm_start_bound_push' : 1e-9,
    #         # 'ipopt.warm_start_bound_frac' : 1e-9,
    #         # 'ipopt.warm_start_slack_bound_frac' : 1e-9,
    #         # 'ipopt.warm_start_slack_bound_push' : 1e-9,
    #         # 'ipopt.warm_start_mult_bound_push' : 1e-9,
    #         'ipopt.dual_inf_tol': 2.,
    #         'ipopt.constr_viol_tol' : 5e-1,
    #         'ipopt.compl_inf_tol' : 5e-1,
    #         'ipopt.acceptable_tol' : 1.5e-1,
    #         # 'ipopt.acceptable_constr_viol_tol' :0.02,
    #         # 'ipopt.acceptable_dual_inf_tol ' : 1e10,
    #         # 'ipopt.acceptable_compl_inf_tol ' : 0.02,
    #         'ipopt.mu_strategy' : 'adaptive'}
    
    opts = {'ipopt.print_level': 1,
            'ipopt.tol' : 4e-3, 
            'ipopt.max_cpu_time': 0.1,
            'ipopt.warm_start_init_point' : 'yes',
            'ipopt.dual_inf_tol': 2.,
            'ipopt.constr_viol_tol' : 5e-1,
            'ipopt.compl_inf_tol' : 5e-1,
            'ipopt.acceptable_tol' : 1.5e-1,
            # 'ipopt.mu_strategy' : 'adaptive'
            }
    solver = nlpsol('solver', 'ipopt', prob, opts)


    w0 = np.concatenate(w0, axis = None)
    lbw = np.concatenate(lbw, axis = None)
    ubw = np.concatenate(ubw, axis = None)
    lbg = np.concatenate(lbg, axis = None)
    ubg = np.concatenate(ubg, axis = None)
    # sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    # w_opt = sol['x'].full().flatten()
    return [solver, lbw, ubw, lbg, ubg]




def fit_poly_44(x,y,z):
    A_matrix = np.vstack([np.ones([x.shape[0]]), x, y, pow(x, 2), x*y, pow(y, 2), pow(x, 3), pow(x,2)*y, pow(y,2)*x, pow(y,3), pow(x, 4), pow(x, 3)*y, pow(x, 2)*pow(y, 2), x*pow(y,3), pow(y, 4)])
    # A_matrix = np.vstack([np.ones([x.shape[0]]), x, y, pow(x, 2), x*y, pow(y, 2)])
    A_matrix = A_matrix.T
    s_vec = np.matrix(z).T
    p = np.linalg.pinv(A_matrix)@s_vec
    return p

# def fit_poly_33(x,y,z):
#     A_matrix = np.vstack([np.ones([x.shape[0]]), x, y, pow(x, 2), x*y, pow(y, 2), pow(x, 3), pow(x,2)*y, pow(y,2)*x, pow(y,3) ])
#     # A_matrix = np.vstack([np.ones([x.shape[0]]), x, y, pow(x, 2), x*y, pow(y, 2)])
#     A_matrix = A_matrix.T
#     s_vec = np.matrix(z).T
#     p = np.linalg.pinv(A_matrix)@s_vec

#     # for i in range(p.shape[0]):
#     #     if i == 0:
#     #         result_p = np.array(p[i])
#     #     else:
#     #         result_p = np.hstack([result_p, p[i]])

#     # print(result_p)
#     return p



def getTerParams(w_opt):
    x_opt = w_opt[1::9]
    y_opt = w_opt[2::9]
    psi_opt = w_opt[3::9]
    v_opt = w_opt[4::9]
    r_opt = w_opt[5::9]
    ux_opt=  w_opt[6::9]
    sa_opt = w_opt[7::9]
    sr_opt = w_opt[8::9]
    ax_opt = w_opt[9::9]
    # params = np.empty((5), float)
    half_length = 10.0
    sampleNum = 20
    x_sample = np.linspace(-half_length, half_length, sampleNum)
    y_sample = np.linspace(-half_length, half_length, sampleNum)
    fitting_num = sampleNum * sampleNum
    for i in range(10+1):
        x = x_opt[i]
        y = y_opt[i]

        lane_length_fitting_x = np.zeros(fitting_num)
        lane_length_fitting_y = np.zeros(fitting_num)
        lane_length_fitting_z = np.zeros(fitting_num)

        count = 0
        for x_idx in range(sampleNum):
            for y_idx in range(sampleNum):
                new_x = x + x_sample[x_idx]
                new_y = y + y_sample[y_idx]
                new_z = -9*np.sin(new_x/15)*np.sin(new_y/15)
                lane_length_fitting_x[count] = new_x
                lane_length_fitting_y[count] = new_y
                lane_length_fitting_z[count] = new_z

        param_fitted = fit_poly_44(lane_length_fitting_x,lane_length_fitting_y,lane_length_fitting_z)
        param_fitted = np.array(param_fitted.T)[0]
        # print(np.array(param_fitted.T)[0])

        # zparx = -0.4 * cos(x/15) * sin(y/15)*9/6
        # zpary = -0.4 * sin(x/15) * cos(y/15)*9/6
        # zparxx = 0.4/15*sin(x/15)*sin(y/15)*9/6
        # zparyy = 0.4/15*sin(x/15)*sin(y/15)*9/6
        # zparxy = -0.4/15*cos(x/15)*cos(y/15)*9/6

        if i == 0:
            params = param_fitted
        else:
            params = np.hstack((params, param_fitted))

    # print(params)
    return params


def main():
    # Terrain Parameter:
    # Specifty OCP parameters
    nu = 2
    nx = 7
    npar = 5
    OCPParams = {
        "XL":[-2, -2, -pi/2, -2, -pi/2, 1, -pi/9],
        "XU":[102, 102, pi/2, 2, pi/2, 9, pi/9],   
        "CL":[-0.5, -2.5],
        "CU":[0.5, 2.5],
        "XF_S":3.6,
        "OBS_SMX":1.0,
        "OBS_SMY":1.0, 
        "nx":nx,
        "nu":nu,
        "npar":npar,
        "T":3,
        "Nck":10,
        "IntegrationScheme":"Trapezoidal"
    }

    zparx = 0.0
    zpary = 0.3

    # zpary = 0.0
    zparxx = 0.0
    zparyy = 0.0
    zparxy = 0.0

    T =  OCPParams["T"]
    goal = [100, 100]
    obslist = np.matrix([[50, 50, 25, 25],[15,15,8,8],[70, 10, 10, 10],[10, 70, 10, 10],[50, 90, 10, 10],[90, 50, 10, 10]])
    init_states = np.array([0.1, 0.1 , 0.01, 0.1, 0.1, 5.1, 0.01])

    w_opt = np.hstack([1e-3,  np.tile(np.hstack([init_states, np.array([0.01, 0.01])]), OCPParams["Nck"]+1)])  #np.zeros(Nck*(7+2)+1)
    w_opt[0] = 0.0


    params = getTerParams(w_opt)
    # params = np.tile(np.array([zparx, zpary, zparxx, zparyy,zparxy]), OCPParams["Nck"]+1)
    
    solver, lbw, ubw, lbg, ubg = defineSolver(init_states, goal, obslist, OCPParams)
###########################

    
    dt = 0.001
    detailed_time = np.arange(0, T, dt)
    data_size = detailed_time.shape[0]
    mat = np.zeros([data_size, 3])

    w_opt_list = np.empty((0,100), float)
    solve_time_list = np.empty((0, 1), float)

    

    while(True):
        if os.path.exists("stop.txt"):
            os.remove("stop.txt")
            break
        while(True):
            if os.path.exists("stop.txt"):
                break
            if os.path.exists("MRZR_Controller_cur_states.txt"):
                sim_cur_states = loadtxt("MRZR_Controller_cur_states.txt", delimiter="\t", unpack=False)
                init_states = sim_cur_states[1:]
                os.remove("MRZR_Controller_cur_states.txt")
                time.sleep(0.001)


                # print(mat)
                np.savetxt('E:/workspace/Vehicle_Modeling_CPP/MRZR_Control/bin/Release/MPC_cmd.txt',mat ,fmt='%.3f')
                time.sleep(0.001)
                np.savetxt('E:/workspace/Vehicle_Modeling_CPP/MRZR_Control/bin/Release/cmd_sent.txt', [0] ,fmt='%.3f')
                time.sleep(0.001)
                predicted_init_states = init_states_predict(init_states, mat[:100,:])

                ocp_start = time.time()
                # w_opt = runMPC(predicted_init_states, goal, obslist, ter_params, f, OCPParams)

                nx = OCPParams["nx"]
                lbw[1:nx+1] = ubw[1:nx+1] = predicted_init_states # current iteration
                # w_opt = np.hstack([1e-3,  np.tile(np.hstack([predicted_init_states.T, np.array([1e-3, 0])]), OCPParams["Nck"]+1)])  #np.zeros(Nck*(7+2)+1)


                # params = np.tile(np.array([zparx, 0.4, zparxx, zparyy,zparxy]), OCPParams["Nck"]+1)
                w_opt[1:nx+1] = predicted_init_states
                params = getTerParams(w_opt)


                sol = solver(lbx = lbw,
                        ubx = ubw,
                        x0  = w_opt,
                        p = params,
                        lbg = lbg,
                        ubg = ubg)
                w_opt = sol["x"].full().flatten()
                # params = getTerParams(w_opt)


                ocp_end = time.time()
                solve_time_list = np.vstack((solve_time_list, ocp_end-ocp_start))
                w_opt_list = np.vstack((w_opt_list, w_opt))

                plt.figure(1)
                plt.clf()
                plt.plot(w_opt[1::9], w_opt[6::9], marker='o', color = 'Red')
                plt.grid(True)
                plt.xlabel('time [s]')
                plt.ylabel('sa [m/s]')
                plt.title('optimal steering angle')
                plt.pause(0.001)
                print("ocp time:")
                print(-ocp_start+ocp_end)
                mat = TransferData(w_opt, OCPParams)
                time.sleep(0.001)
                break
            else:
                time.sleep(0.01)

    print(np.max(solve_time_list))
    print(np.mean(solve_time_list))
    
    KeySolutions = {
        "Solution":w_opt_list,
        "Time":solve_time_list
    }

    from scipy.io import savemat
    savemat("Result_Mat\OneObsChallengeTerrain_Ter.mat", KeySolutions)

main()