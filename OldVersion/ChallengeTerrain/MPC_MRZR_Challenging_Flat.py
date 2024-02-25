from casadi import *
import matplotlib.pyplot as plt
import copy
from numpy import loadtxt
import os.path
import time


def mapfunc(X, Y):
    return -6*np.sin(X/15)*np.sin(Y/15)+6

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



    normal_vec  = np.array([0, 0, 1])
    normal_vec  = normal_vec/np.linalg.norm(normal_vec)
    theta       = np.arcsin(normal_vec[0])
    phi         = -np.arcsin(normal_vec[1]/np.arccos(theta))
    RX          = np.array([[1.0,0.0,0.0], [0.0, np.cos(phi), -np.sin(phi)],[0.0, np.sin(phi), np.cos(phi)]])
    RY          = np.array([[np.cos(theta), 0, np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)] ]  )  
    RZ          = np.array([[np.cos(psi),-np.sin(psi), 0],[np.sin(psi),np.cos(psi),0],[0, 0, 1] ]  )   
    R           = RX@ RY@RZ 
    ERX         = np.array([[1,0,np.sin(theta)],[0,np.cos(phi), -np.cos(theta)*np.sin(phi)],[0,np.sin(phi),np.cos(phi)*np.cos(theta)]])



    wx          = 0
    wy          = 0

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
    KFZYF       = 318.2923 
    KFZYR       = 406.5813 
    rho_ks      = 2.5 
    TS_FZ_circle= 100
    TS_FZ       = 400
    ax_max      = 3.0
    ax_min      = -3.0
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



    normal_vec  = horzcat(0, 0 , 1)
    normal_vec  = normal_vec/norm_1(normal_vec)
    theta       = asin(normal_vec[0])
    phi         = -asin(normal_vec[1]/acos(theta))

    RX          = vertcat(horzcat(1.0,0.0,0.0), horzcat(0.0, cos(phi), -sin(phi)),horzcat(0.0, sin(phi), cos(phi)))
    RY          = vertcat( horzcat(cos(theta), 0, sin(theta)), horzcat(0,1,0), horzcat(-sin(theta),0,cos(theta))  )  
    RZ          = vertcat( horzcat(cos(psi),-sin(psi), 0), horzcat(sin(psi),cos(psi),0), horzcat(0, 0, 1) )   
    R           = RX@ RY @RZ 
    ERX         = vertcat( horzcat(1,0,sin(theta)), horzcat(0,cos(phi), -cos(theta)*sin(phi)),horzcat(0,sin(phi),cos(phi)*cos(theta)))

    wx          = 0
    wy          = 0

    G           = R.T @ vertcat(0,0,-g)
    gx          = G[0]
    gy          = G[1]
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
    dpsi        = r
    dv          = (FY1 + FY2) / M -  r * ux + gy
    dr          = (FY1 * la - FY2 * lb) / Izz
    du          = ax
    dsa         = sr
    dstates     = vertcat(dx, dy, dpsi, dv, dr, du, dsa)

    FZFL_con    = FZF / 2 - ((FY1 + FY2) / M) * KFZYF - TS_FZ
    FZFR_con    = FZF / 2 + ((FY1 + FY2) / M) * KFZYF - TS_FZ
    FZRL_con    = FZR / 2 - ((FY1 + FY2) / M) * KFZYR - TS_FZ
    FZRR_con    = FZR / 2 + ((FY1 + FY2) / M) * KFZYR - TS_FZ
    # FORCE_con   =  ((ax - r*v - gx) * M )**2 -  ( ( 2 * KFZF * gz + 2 * KFZR * gz) * mu)**2 + TS_FZ_circle**2 
    # FORCE_con   = 0
    ax_max_con  = ax - r*v - gx - ax_max 
    ax_min_con  = ax_min - ax + r * v + gx
    # ks_con         = exp(rho_ks * (- FZFL_con)) + exp(rho_ks * (- FZFR_con)) + exp(rho_ks * (- FZRL_con)) + exp(rho_ks * (- FZRR_con))+ exp(rho_ks * ax_max_con) + exp(rho_ks * ax_min_con)
    ks_con      = 0
    con = vertcat(ks_con, 0)
    L = 0.1 * sa**2 + 0.05 * ax**2 + 0.12 * v**2 + 0.3 * sr**2 + 0.02 * (alpha_f**2 + alpha_r**2) + 0.1 * gy**2
    return dstates, con, L

def runMPC(init_states, goal, obslist, ter_params, f, OCPParams):
    goal_tolerance = OCPParams["XF_S"]
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
    mystates = init_states


    lbw.append(mystates)
    ubw.append(mystates)
    w0.append(mystates)
    Uk = MX.sym('U_0', nu)
    w.append(Uk)
    lbw.append(CL)
    ubw.append(CU)
    w0.append([0, 0])



    [fj, conj, qj] = f(Xk, Uk, ter_params)

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
        [fjpl, conjpl, qjpl] = f(Xkpl, Ukpl, ter_params)
        J = J + qjpl*dt


        if (OCPParams["IntegrationScheme"] == "Trapezoidal"):
            g.append(Xkpl - (fjpl+fj)/2 *dt -Xk)
        elif (OCPParams["IntegrationScheme"] == "BkwEuler"):
            g.append(Xkpl - (fjpl)/2 *dt -Xk)
        elif (OCPParams["IntegrationScheme"] == "MidCol"):
            [fjmid, conjmid, qjmid] = f(1/2*(Xk+Xkpl)+dt/8*(fj-fjpl), 1/2*(Uk+Ukpl), ter_params)
            g.append(-3/2*(Xk-Xkpl)-dt/4*(fj+fjpl)-dt*fjmid)
        else:
            raise Exception("No integration scheme available")

        lbg.append(np.zeros(nx))
        ubg.append(np.zeros(nx))

        g.append(conj)
        lbg.append([0.0, -inf])
        ubg.append([1.0, 0.0])

        for obs_index in range(obslist_sm.shape[0]):
            g.append((   (Xkpl[0]-obslist_sm[obs_index,0])/obslist_sm[obs_index,2]   )**2+( (Xkpl[1]-obslist_sm[obs_index,1])/obslist_sm[obs_index,3] )**2)
            lbg.append(1.0)
            ubg.append(inf)

        Xk = Xkpl
        Uk = Ukpl
        fj = fjpl
        qj = qjpl
        conj = conjpl


    J = J+(1e3)*mu**2 + 150.0*(  ((Xk[0]-goal[0])**2+(Xk[1] - goal[1])**2)/((mystates[0]-goal[0])**2+(mystates[1] - goal[1])**2)    )
    w = vertcat(*w)
    g = vertcat(*g)
    prob = {'f': J, 'x': w, 'g': g}
    opts = {'ipopt.print_level': 0}
    solver = nlpsol('solver', 'ipopt', prob, opts)
    w0 = np.concatenate(w0, axis = None)
    lbw = np.concatenate(lbw, axis = None)
    ubw = np.concatenate(ubw, axis = None)
    lbg = np.concatenate(lbg, axis = None)
    ubg = np.concatenate(ubg, axis = None)
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    w_opt = sol['x'].full().flatten()

    return w_opt

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
    data_size = 500
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
    return states_list[499,:]


def main():
    # Terrain Parameter:
    p00 = 6
    p10 = -3.7328e-06
    p01 = 3.5151e-06
    p20 = -2.7239e-06
    p11 = -0.0266
    p02 = -2.6646e-06
    p30 = 2.5426e-08
    p21 =  -8.4099e-08
    p12 = 1.0947e-07
    p03 = -2.9803e-08
    p40 =  2.2629e-08
    p31 =  1.8947e-05
    p22 = 1.2006e-07
    p13 =  1.8948e-05
    p04 =  2.1190e-08

    ter_params = [p00, p10, p01, p20, p11, p02, p30, p21, p12, p03, p40, p31, p22, p13, p04] 

    # Declare model variables
    x1 = SX.sym('x1')
    x2 = SX.sym('x2')
    x3 = SX.sym('x3')
    x4 = SX.sym('x4')
    x5 = SX.sym('x5')
    x6 = SX.sym('x6')
    x7 = SX.sym('x7')
    x = vertcat(x1, x2, x3, x4, x5, x6, x7)
    nu = 2
    u1 = SX.sym('u1')
    u2 = SX.sym('u2')
    u = vertcat(u1, u2)
    param = SX.sym('ter', 15)
    xdot, con, L_0 = VehicleDynamics3D(x, u, param)
    nx = 7

    L = L_0
    f = Function('f', [x, u, param], [xdot, con, L], ['x', 'u', 'ter'], ['xdot', 'con','L'])

    # Specifty OCP parameters
    OCPParams = {
        "XL":[-105, -105, -pi/2, -2, -pi/2, 1, -pi/9],
        "XU":[105, 105, pi/2, 2, pi/2, 9, pi/9],   
        "CL":[-0.5, -2.5],
        "CU":[0.5, 2.5],
        "XF_S":3.6,
        "OBS_SMX":1.0,
        "OBS_SMY":1.0, 
        "nx":nx,
        "nu":nu,
        "T":3,
        "Nck":10,
        "IntegrationScheme":"Trapezoidal"
    }
    #  Specify goal position
    # goal = np.array([100.0, 100.0])
    
    goal_list = np.array([[100, 100]])
    # goal_list = np.array([[18, 40], [50, 80], [100, 100]])


    goal = goal_list[0, :]
    goal_check = 1
    #  Specify Obstacle Position
    obslist = np.matrix([[50, 50, 25, 25],[15,15,8,8],[70, 10, 10, 10],[10, 70, 10, 10],[50, 90, 10, 10],[90, 50, 10, 10]])


    T =  OCPParams["T"]
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



                np.savetxt('E:/workspace/Vehicle_Modeling_CPP/MRZR_Control/bin/Release/MPC_cmd.txt',mat ,fmt='%.3f')
                time.sleep(0.001)
                np.savetxt('E:/workspace/Vehicle_Modeling_CPP/MRZR_Control/bin/Release/cmd_sent.txt', [0] ,fmt='%.3f')
                time.sleep(0.001)
                predicted_init_states = init_states_predict(init_states, mat[:500,:])

                ocp_start = time.time()
                w_opt = runMPC(predicted_init_states, goal, obslist, ter_params, f, OCPParams)
                ocp_end = time.time()
                solve_time_list = np.vstack((solve_time_list, ocp_end-ocp_start))
                w_opt_list = np.vstack((w_opt_list, w_opt))
                if (((goal[0] - predicted_init_states[0])**2 + (goal[1] - predicted_init_states[1])**2) <= 15**2):
                    # goal_check = goal_check + 1
                    if (goal_check + 1 <= goal_list.shape[0]):
                        goal_check = goal_check+1
                    goal = goal_list[goal_check - 1, :]

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

    print(w_opt_list)
    print(solve_time_list)
    KeySolutions = {
        "Solution":w_opt_list,
        "Time":solve_time_list
    }

    from scipy.io import savemat
    savemat("Result_Mat\ChallengeTerrain_Plane.mat", KeySolutions)
main()