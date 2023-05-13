import numpy as np

# Clean Configuration

rho = .0000491386 #slug/ft^3 @28,000ft
V = 464.148 #ft/s @cruise
u = 0.3251*10**-6 #slug/ft·s @28,000ft
M = 0.423 # @cruise


def skinfriction_drag(rho, V, L, u, M, plff, plfw, struct):
    # rho = density
    # V = velocity
    # L = length of fuselage or chord length
    # u = viscosity

    R = (rho * V * L) / u
    R_cutoff = 38.21*(L/0.0000208)**1.053

    # M = mach number
    # plff = predicted percent of laminar flow over fuselage
    # plff = predicted percent of laminar flow over wing
    # struct = "wing" or "fuselage"

    C_fc_lam = 1.328 / np.sqrt(R)
    if R_cutoff < R:
        C_fc_turb = 0.455 / (((np.log10(R_cutoff))**2.58) * ((1 + 0.144 * (M**2))**0.65))
    else:
        C_fc_turb = 0.455 / (((np.log10(R))**2.58) * ((1 + 0.144 * (M**2))**0.65))

    if struct == "wing":
        Cfc_weightedavg = ((C_fc_lam * plfw)+ (C_fc_turb * (1 - plfw))) / 2
    elif struct == "fuse":
        Cfc_weightedavg = ((C_fc_lam * plff)+ (C_fc_turb * (1 - plff))) / 2

    return Cfc_weightedavg

def wing_formfactor_drag(xc_max, t, c, M, maxsweep):
    # xc_max = location of maximum thickness on wing, normalized by chord length
    # t = maximum thickness
    # c = chord length
    # M = mach number
    # maxsweep = sweep angle of the maximmnum thickness line (should be about the same as the sweep angle if the same airfoil is used all the way)

    FFc = ((1 + (0.6 / xc_max) * (t/c)) + 100*((t/c)**4)) * ((1.34 * (M ** 0.18)) * ((np.cos(maxsweep))**0.28))
    return FFc

def finenessratio(l, d_Amax, cyl):
    # L = length of object
    # d_Amax = either the diameter (for cylindrical / spherical objects) or the maximum cross-sectional area (for noncylindrical objects)
    # cyl = 1, if cylindrical, = 0  if not

    if cyl == 1:
        f = l / d_Amax
    elif cyl !=1:
        f = l / np.sqrt((4/np.pi) * d_Amax)
    return f

def fusenacelle_formfactor_drag(f, main):
    # f = fineness ratio
    # main = 1, if fuselage, = 0 if external (external storage, nacelles, etc)

    if main == 1:
        if f < 6:
            FFc = (0.9 + (5/(f**1.5)) + (f/400))
        elif f >= 6:
            FFc = (1 + (60/(f**3)) + (f/400))
    elif main != 1:
        FFc = 1 + (0.35/f)
    return FFc

def upsweep_drag(ua, A_max, S_ref):
    # ua = upsweep angle
    # A_max = cross-section maximum area
    # S_ref = reference area

    D_q = 3.83 * (ua**2.5) * A_max
    CD_misc = D_q / S_ref
    return CD_misc

def base_drag(M, A_base, S_ref):
    # M = mach number
    # A_base = flat area
    # S_ref = reference area

    D_q = (0.139 + 0.419 * ((M-0.161)**2)) * A_base
    CD_misc = D_q / S_ref
    return CD_misc

def prop_drag(A_disk, A_blade, S_ref, prop):
    # A_disk = total disk area
    # A_blade = blade area
    # S_ref = reference area
    # prop = "feathered" or "stopped"
    
    sigma = A_blade / A_disk
    if prop == 'feathered':
        D_q = 0.1 * sigma * A_disk
    elif prop == 'stopped':
        D_q = 0.8 * sigma * A_disk
    CD_misc = D_q / S_ref
    return CD_misc



S_ref = 923 #ft^2  from 'wing area'#NEED UDPATE

############################## Fuselage ###############################
Swet_fuse = 2035.836 # wetted area of fuselage # NEED UPDATE

L = 83.0 # length ft
D = 9.0 # diameter ft
plff = 0.10
plfw = 0.35
struct = 'fuse'
Cf_fuse = skinfriction_drag (rho, V, L, u, M, plff, plfw, struct)

cyl = 1
f_fuse = finenessratio(L, D, cyl)

main = 1
FF_fuse = fusenacelle_formfactor_drag(f_fuse, main)

Q_fuse = 1.0

A_max = np.pi * ((D/2)**2)
ua = 10.6794 * (np.pi/180) # upsweep angle in radians
misc_upsweep  = upsweep_drag (ua, A_max, S_ref)

CD_fuse = (Cf_fuse * FF_fuse * Q_fuse * Swet_fuse)



########################## Wing ##########################

Swet_wing = 2 * 934.698 # wetted area of wing #NEED UPDATE

L = 8.92 # characteritic length ft

plff = 0.10
plfw = 0.35
struct = 'wing'
Cf_wing = skinfriction_drag (rho, V, L, u, M, plff, plfw, struct)

xc_max = 0.30  # chord normalized location of max thickness # NEED UPDATE
t = 1.1958 # max thickness # NEED UPDATE
c = 8.92 # chord length ft # NEED UPDATE
maxsweep = 19 # sweep angle of max thickness deg line #NEED UPDATE
FF_wing = wing_formfactor_drag(xc_max, t, c, M, maxsweep)

Q_wing = 1.0

CD_wing = (Cf_wing * FF_wing * Q_wing * Swet_wing)


################################ Nacelles ###################################
Swet_nacelle = 2 * 315.23 # wetted area of nacelle ft^2 # NEED UPDATE

L = 14.8 # length ft
D = 2.1142 # diameter ft
plff = 0.10
plfw = 0.35
struct = 'fuse'
Cf_nacelle = skinfriction_drag (rho, V, L, u, M, plff, plfw, struct)

cyl = 1
f_nacelle = finenessratio(L, D, cyl)

main = 0
FF_nacelle = fusenacelle_formfactor_drag(f_nacelle, main)

Q_nacelle = 1.5

d_prop = 12 #ft
A_disk = np.pi * (d_prop/2) ** 2 #ft^2
A_blade = (A_disk / 22) #NEED UPDATE
prop = 'feathered'
misc_prop = prop_drag(A_disk, A_blade, S_ref, prop)

CD_nacelle = 2 * (Cf_nacelle * FF_nacelle * Q_nacelle * Swet_nacelle)




############################### Vertical Tail #################################
Swet_vtail = 193 # wetted area of vertical tail # NEED UPDATE

c = 14.77 # chord length # NEED UPDATE

L = c # characteritic length ft

plff = 0.10
plfw = 0.35
struct = 'wing'
Cf_vtail = skinfriction_drag (rho, V, L, u, M, plff, plfw, struct)

xc_max = 0.30 # chord normalized location of max thickness # NEED UPDATE
t = 1.333 # max thickness # NEED UPDATE
maxsweep = 25.3 # sweep angle of max thickness line deg # NEED UPDATE
FF_vtail = wing_formfactor_drag(xc_max, t, c, M, maxsweep)

Q_vtail = 1.05

CD_vtail = (Cf_vtail * FF_vtail * Q_vtail * Swet_vtail)



############################### Horizontal Tail #################################
Swet_htail = 2 * 218 # wetted area of horizontal tail # NEED UPDATE

c = 7.38 # chord length # NEED UPDATE

L = c # characteritic length ft

plff = 0.10
plfw = 0.35
struct = 'wing'
Cf_htail = skinfriction_drag (rho, V, L, u, M, plff, plfw, struct)

xc_max = 0.30 # chord normalized location of max thickness # NEED UPDATE
t = 0.8856 # max thickness # NEED UPDATE
maxsweep = 30 # sweep angle of max thickness line
FF_htail = wing_formfactor_drag(xc_max, t, c, M, maxsweep)

Q_htail = 1.05

CD_htail = (Cf_htail * FF_htail * Q_htail * Swet_htail)


################## Miscellaneous ###################

CD_misc = misc_upsweep + misc_prop 

A_front_wheel = 2*2.875*0.85 # NEEDS UPDATE
CD_wheel = (0.25 * A_front_wheel) / S_ref

A_front_wheel2 = 4*3.6*1.1 # NEEDS UPDATE
CD_wheel2 = (0.25 * A_front_wheel2) / S_ref

A_front_strut =   (2.350+1.5+2.4)*0.3 # NEEDS UPDATE
CD_strut = (0.30 * A_front_strut) / S_ref

A_front_Tflaps = 140.666 * 0.4 # ft^2 # NEEDS UPDATE
CD_Tflaps = (1.28 * A_front_Tflaps) / S_ref

A_front_Lflaps = 140.666 * 0.4 # ft^2 # NEEDS UPDATE
CD_Lflaps = (1.28 * A_front_Lflaps) / S_ref

######################### Flight Configurations ################################

# Clean Configuration
CD_0_Clean = (((CD_fuse + CD_wing + CD_nacelle + CD_vtail + CD_htail) / S_ref) + CD_misc) * 1.05 #extra 5% accounts for leakages and protuberances
print(CD_0_Clean)


# Takeoff Flaps, Landing Gear Up
CD_0_TFGU = (((CD_fuse + CD_wing + CD_nacelle + CD_vtail + CD_htail) / S_ref) + CD_misc + CD_Tflaps) * 1.05
print("Takeoff flaps, LG up " + str(CD_0_TFGU))

# Takeoff Flaps, Landing Gear Down
CD_0_TFGD = (((CD_fuse + CD_wing + CD_nacelle + CD_vtail + CD_htail) / S_ref) + CD_misc + CD_Tflaps + CD_wheel + CD_wheel2 + CD_strut) * 1.05
print("Takeoff flaps, LG down " + str(CD_0_TFGD))


# Landing Flaps, Landing Gear Up
CD_0_LFGU = (((CD_fuse + CD_wing + CD_nacelle + CD_vtail + CD_htail) / S_ref) + CD_misc + CD_Lflaps) * 1.05
print("ladning flaps, LG up " + str(CD_0_LFGU))

# Landing Flaps, Landing Gear Down
CD_0_LFGD = (((CD_fuse + CD_wing + CD_nacelle + CD_vtail + CD_htail) / S_ref) + CD_misc + CD_Lflaps + CD_wheel + CD_wheel2 + CD_strut) * 1.05
print("ladning flaps, LG down " + str(CD_0_LFGD))