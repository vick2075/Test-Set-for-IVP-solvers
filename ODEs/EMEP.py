import numpy as np

# EMEP Problem
        
tt = np.empty(10)
tt[0] = 4*3600
tt[1] = 20*3600
for i in range(1,5):
    tt[2*i] = tt[0] +24*3600*i
    tt[2*i+1] = tt[1] +24*3600*i
    

# listing of species

ls = ['NO',     'NO2',    'SO2',    'CO',     'CH4',     'C2H6',
     'NC4H10', 'C2H4',   'C3H6',   'OXYL',   'HCHO',    'CH3CHO',
     'MEK',    'O3',     'HO2',    'HNO3',   'H2O2',    'H2',
     'CH3O2',  'C2H5OH', 'SA',     'CH3O2H', 'C2H5O2',  'CH3COO',
     'PAN',    'SECC4H', 'MEKO2',  'R2OOH',  'ETRO2',   'MGLYOX',
     'PRRO2',  'GLYOX',  'OXYO2',  'MAL',    'MALO2',   'OP',
     'OH',     'OD',     'NO3',    'N2O5',   'ISOPRE',  'NITRAT',
     'ISRO2',  'MVK',    'MVKO2',  'CH3OH',  'RCO3H',   'OXYO2H',
     'BURO2H', 'ETRO2H', 'PRRO2H', 'MEKO2H', 'MALO2H',  'MACR',
     'ISNI',   'ISRO2H', 'MARO2',  'MAPAN',  'CH2CCH3', 'ISONO3',
     'ISNIR',  'MVKO2H', 'CH2CHR', 'ISNO3H', 'ISNIRH',  'MARO2H']



# mixing height in cm
HMIX = 1.2e+5

# distribution of VOC emissions among species
FRAC6=0.07689
FRAC7=0.41444
FRAC8=0.03642
FRAC9=0.03827
FRAC10=0.24537
FRAC13=0.13957

VEKT6=30.0
VEKT7=58.0
VEKT8=28.0
VEKT9=42.0
VEKT10=106.0
VEKT13=46.0

VMHC=1.0/(FRAC6/VEKT6+FRAC7/VEKT7+FRAC8/VEKT8+ \
          FRAC9/VEKT9+FRAC10/VEKT10+FRAC13/VEKT13)


# values for NOX and HC emissions in molecules/cm**2xs
EMNOX = 2.5e+11
EMHC = 2.5e+11


# rural case
FACISO = 1.0
FACHC = 1.0
FACNOX = 1.0

# Emissions in molec/(cm^2*s) of NO, NO2, SO2, CO, CH4, C2H6, NC4H10,
#                                C2H4, C3H6, O-XYLENE, C2H5OH and ISOPRENE
EMIS1 = EMNOX * FACNOX
EMIS2 = 0.0
EMIS3 = EMNOX * FACNOX
EMIS4 = EMHC*10.0 * FACHC
EMIS5 = 0.0
EMIS6 = EMHC*FRAC6/VEKT6*VMHC  * FACHC
EMIS7 = EMHC*FRAC7/VEKT7*VMHC  * FACHC
EMIS8 = EMHC*FRAC8/VEKT8*VMHC  * FACHC
EMIS9 = EMHC*FRAC9/VEKT9*VMHC  * FACHC
EMIS10= EMHC*FRAC10/VEKT10*VMHC* FACHC
EMIS11= 0.5 * FACISO * EMHC
EMIS12= 0.0
EMIS13= EMHC*FRAC13/VEKT13*VMHC* FACHC


# time-dependent EMEP coefficients
M = 2.55e+19
O2 = 5.2e+18
XN2 = 1.99e+19
RPATH3 = 0.65
RPATH4 = 0.35


# ODEs

def f(t, u):
    # dissociation rate coefficient
    TIME = t # to input
    TIMEH=TIME/3600.
    ITIMEH=int(TIMEH)
    I24HRS=ITIMEH//24+1
    TIMEOD= TIMEH-(I24HRS-1)*24.0 # int(24. -(I24HRS-1))+(TIMEH-ITIMEH)#

    # XLHA local hour angle
    XLHA=(1.0+TIMEOD*3600.0/4.32e+4)*np.pi

    # FI (Norwegian PHI!) latitude, dekl solar declination
    # here latitude = 50 deg. N
    FI=50.0*np.pi/180.0
    DEKL=23.5*np.pi/180.0
    SEC=1.0/(np.cos(XLHA)*np.cos(FI)*np.cos(DEKL)+np.sin(FI)*np.sin(DEKL))


    #  def of temperature variation
    #     XP=8.7D-5*TIMEOD*3600.D0-2.83
    #     T=8.3D0*SIN(XP)+289.86
    #   for simplicity
    # T temperature in K
    T=298.0

    #   def of water vapor concentration
    XQ=-7.93e-5*TIMEOD*3600.0+2.43
    RH=23.0*np.sin(XQ)+66.5

    XZ=(597.3-0.57*(T-273.16))*18.0/1.986*(1.0/T-1.0/273.16)
    H2O=6.1078*np.exp(-XZ)*10.0/(1.38e-16*T)*RH




    A = [1.23e-3,2.00e-4,1.45e-2,2.20e-5,3.00e-6, \
         5.40e-5,6.65e-5,1.35e-5,2.43e-5,5.40e-4, \
         2.16e-4,5.40e-5,3.53e-2,8.94e-2,3.32e-5,2.27e-5]

    B = [0.60,   1.40,   0.40,   0.75,   1.25, \
         0.79,   0.60,   0.94,   0.88,   0.79, \
         0.79,   0.79,  0.081,  0.059,   0.57,   0.62]



    # SEC = 1./np.cos(THETA) # THETA is solar zenith angle


    # Calculate  values of photolysis rates DJ(1..16), based
    #   upon RGD A & B coefficients and correction factors from HOUGH (1988)
    DJ = np.empty(16)

    if TIMEOD < 4.0 or TIMEOD >= 20.0:
        # in the dark
        for i in range(16):
            DJ[i] = 1.0e-40
    else:
        # daytime
        for i in range(16):
            DJ[i] = A[i]*np.exp(-B[i]*SEC)
            if DJ[i] < 0.0:
                break


    # Set up chemical reaction rate coefficients:

    RC = np.empty(266)

    for i in range(266):
        RC[i] = 0.0

    # assuming 80% N2, 20% O2
    RC[0]=5.7e-34*(T/300.0)**(-2.8) * M * O2

    # unchanged, as source of O+NO+O2 reaction rate unknown.
    RC[4]=9.6e-32*(T/300.0)**(-1.6) * M

    # estimate from Demore(1990) (=A92) O(d)+N2 and DeMore O(D)+O2
    RC[6]  =2.0e-11*np.exp(100.0/T) * M

    RC[10] =1.8e-12*np.exp(-1370.0/T)
    RC[11] =1.2e-13*np.exp(-2450.0/T)
    RC[12] =1.9e-12*np.exp(-1000.0/T)
    RC[13] =1.4e-14*np.exp(-600.0/T)
    RC[14] =1.8e-11*np.exp(+110.0/T)
    RC[16] =3.7e-12*np.exp(240.0/T)

    # M.J (from Wayne et al.)
    RC[18] =7.2e-14*np.exp(-1414.0/T)
    # RC[18] =7.2e-13*np.exp(-1414.0/T)

    # M.J suggests that rc(27) be omitted:
    #     RC[26] =8.5e-13*np.exp(-2450.0/T)
    #  My change to get similar to A92 troe.
    RC[28] =7.1e+14*np.exp(-11080.0/T)

    RC[29] =4.8e-11*np.exp(+250.0/T)

    # De More,1990 .. no change : oh + h2o2
    RC[30] =2.9e-12*np.exp(-160.0/T)

    # : oh + h2
    RC[32] =7.7e-12*np.exp(-2100.0/T)

    # My, similar to DeMore et al complex : oh+hno3
    RC[34] =1.0e-14*np.exp(785.0/T)

    # Mike Jenkin`s suggestion for ho2 + ho2 reactions: (from DeMore et al.)
    RC[35] =2.3e-13*np.exp(600.0/T)
    RC[35] = RC[35] + M * 1.7e-33*np.exp(1000.0/T)
    RC[35] = RC[35] * (1.0 + 1.4e-21 * H2O *np.exp(2200.0/T))

    RC[58] =3.9e-12*np.exp(-1885.0/T)
    # : ch3o2 + no
    RC[59] =4.2e-12*np.exp(180.0/T)
    # A92 + A90 assumption that ka = 0.5D0 * k
    RC[60] =5.5e-14*np.exp(365.0/T)

    # A92 + A90 assumption that kb = 0.5D0 * k
    RC[61] =5.5e-14*np.exp(365.0/T)

    RC[62] =3.3e-12*np.exp(-380.0/T)
    # : ho2 + ch3o2
    RC[64] =3.8e-13*np.exp(780.0/T)

    # new: ch3ooh + oh -> hcho + oh
    RC[66] =1.0e-12*np.exp(190.0/T)
    # new: ch3ooh + oh -> ch3o2
    RC[67] =1.9e-12*np.exp(190.0/T)

    RC[70] =7.8e-12*np.exp(-1020.0/T)
    # new: c2h5o2 + ho2 -> c2h5ooh (r2ooh)
    RC[73] =6.5e-13*np.exp(650.0/T)

    RC[74] =5.6e-12*np.exp(310.0/T)
    # TOR90 assumption w.r.t. rc(67) : c2h5ooh + oh -> ch3cho + oh
    #     RC[75] = 5.8 * RC[66]
    RC[75] = 5.8e-12*np.exp(190.0/T)
    # : approximation to troe expression
    RC[77] =1.34e+16*np.exp(-13330.0/T)
    # additional reactions :-
    # : ho2 + ch3coo2 -> rco3h
    RC[87] =1.3e-13*np.exp(1040.0/T)
    # : ho2 + ch3coo2 -> rco2h + o3
    RC[88] =3.0e-13*np.exp(1040.0/T)

    RC[93] =2.8e-12*np.exp(530.0/T)
    # D & J, gives results very close to A92
    RC[80] =1.64e-11*np.exp(-559.0/T)


    RC[82] =RC[59]
    RC[104]=RC[59]

    RC[109]=RC[59]
    # A90/PS  isnir + no
    RC[161]=RC[59]
    # A90/PS  isono3 + no
    RC[163]=RC[59]

    # From RGD, but very similar to A92 Troe
    RC[108]=1.66e-12*np.exp(474.0/T)

    RC[111]=1.2e-14*np.exp(-2630.0/T)

    RC[122]=6.5e-15*np.exp(-1880.0/T)

    RC[125] = RC[59]

    RC[219] = RC[59]

    RC[235] = RC[59]
    #  natural voc reactions #
    #  isoprene + o3 -> products
    RC[149] = 12.3e-15*np.exp(-2013.0/T)
    #  isoprene + oh -> isro2
    RC[150] = 2.54e-11*np.exp(410.0/T)
    #  isoprene-RO2 + no
    RC[151] = RC[59]
    #  methylvinylketone (mvk) + oh
    RC[152] = 4.13e-12*np.exp(452.0/T)
    #  mvko2 + no
    RC[153] = RC[59]
    #  macr + oh
    RC[157] = 1.86e-11*np.exp(175.0/T)
    #  ch2cch3 + no
    RC[158] = RC[59]
    # mvk + o3
    RC[159] = 4.32e-15*np.exp(-2016.0/T)

    #     RC[254] = 9.9e-16*np.exp(-731.0/T)

    #     aerosol formation and depositions ..x
    #     parameterization of heterogeneous loss in 1./s for humidities
    #     less than 90%.  applied to hno3, h2o2, h2so4, ch3o2h.
    RC[42]=5.e-6
    if(RH > 0.90):
        RC[42]=1.e-4

    RC[43] =RC[42]
    RC[44] =RC[42]


    RC[7]  =2.2e-10
    # My (?) , similar to A92 troe expression.
    RC[19] =1.4e-12
    # My (?) , similar to A92 troe expression.
    RC[20] =1.4e-11
    # RGD Note (=DeMore?)
    RC[25] =4.1e-16
    # RGD
    RC[38] =1.35e-12
    # RGD
    RC[39] =4.0e-17
    # A92, assuming all products are ch3cho
    RC[63] =3.2e-12
    # A92, with temperature dependance neglected.
    RC[65] =9.6e-12

    RC[68] =5.8e-16
    RC[69] =2.4e-13
    RC[71] =8.9e-12

    # A92 : approximation to troe expression
    RC[76] =1.0e-11
    # : ch3coo2 + no
    RC[78] =2.0e-11
    #: sum of two pathways
    RC[79] =1.1e-11
    # cmj   RC[83] =2.5e-14
    #   kho2ro2 : estimate of rate of ho2 + higher RO2, suggested by Yvonne
    RC[84] =1.0e-11
    # ignoring slight temperature dependance
    RC[85] =1.15e-12
    # MJ suggestion.. combine oh + secc4h9o2, meko2 rates   =>
    #     RC[86]= RC[67] + RC[85], approx. at 298
    RC[86]= 4.8e-12
    # new A92 : oh + ch3co2h -> ch3o2
    RC[89] =8.0e-13
    # Approximates to A92 Troe ...
    RC[124]=2.86e-11
    # rate for ch2chr + oh, = k68+k125 (propene), Yv.
    RC[145] = 3.2e-11
    # rate for isono3h + oh, = k156 (isro2h), Yv.
    RC[146] = 2.0e-11
    # MY GUESS rate of oh + mvko2h
    RC[147] = 2.2e-11
    # MY GUESS rate of oh + other biogenic ooh
    RC[148] = 3.7e-11
    #   kho2ro2 : estimate of rate of ho2 + higher RO2, suggested by Yvonne
    #  isro2 + ho2
    RC[154] =RC[84]
    #  isro2h + oh
    RC[155] =2.0e-11
    #   isro2h + o3
    RC[156] =8.0e-18
    #   isni + oh
    RC[160] =3.35e-11
    #  isopre + no3
    RC[162] =7.8e-13
    #  RC[162] =7.8e-16
    # Unchanged, also in IVL scheme
    RC[218]=2.0e-11

    RC[220]=1.1e-11

    RC[221]=1.7e-11
    # MJ suggestion.. combine oh + malo2h rates   =>
    #     RC[222]= RC[67] + RC[218]
    RC[222]= 2.4e-11

    RC[233]=1.37e-11
    # MJ suggestion.. combine rc(68) with rc(234) =>
    #     RC[234]= RC[67] + RC[233]
    RC[234]= 1.7e-11

    #         deposition loss rate coefficients vd/hmix, vd in cm/s.
    #         hno3     calculated     rc(46)
    #         so2      0.5            rc(47)
    #         h2o2     0.5            rc(47)
    #         no2      0.2            rc(48)
    #         o3       0.5            rc(49)
    #         pan      0.2            rc(50)
    #         h2so4    0.1            rc(51)
    # simple approx. for now - reduce all vg by 4 at night to allow
    #    for surface inversion......
    delta=1.0
    if(TIMEOD >= 20.0 or TIMEOD < 4.0):
        delta=0.25
    #         if(timeod.ge.20.D0.or.timeod.le.4.D0) delta=0.25D0
    #         if (night) delta=0.25D0
    RC[45] =2.0 * delta /HMIX
    RC[46] =0.5 * delta /HMIX
    RC[47] =0.2 * delta /HMIX
    RC[48] =0.5 * delta /HMIX
    RC[49] =0.2 * delta /HMIX
    RC[50] =0.1 * delta /HMIX
    # dep. of organic peroxides = 0.5 cms-1
    RC[51] = 0.5 *delta /HMIX
    # dep. of ketones, RCHO  = 0.3 cms-1
    RC[52] = 0.3 *delta /HMIX
    
    s1 = u[37-1]*(RC[13-1]*u[14-1]+RC[31-1]*u[17-1]+RC[33-1]*u[18-1]+RC[39-1]*u[3-1]+RC[63-1]*u[46-1]+RC[64-1]* \
                  u[20-1]+RC[66-1]*u[11-1]+RC[70-1]*u[4-1]+RC[221-1]*u[32-1])+u[19-1]* \
                  (RC[40-1]*u[3-1]+2.0*RC[61-1]*u[19-1]+0.5*RC[80-1]*u[24-1])+DJ[11-1]*u[30-1]+u[1-1]*RC[60-1]* \
                  (u[19-1]+u[29-1]+u[31-1]+u[33-1]+u[35-1]+0.95*u[45-1]+u[26-1]*RPATH3+0.78*u[43-1]+u[59-1]+ \
                   0.5e-1*u[61-1]+0.8*u[60-1])+RC[72-1]*u[1-1]*u[23-1]+DJ[8-1]*u[12-1]
    s2 = 2.0*RC[8-1]*H2O*u[38-1]+u[15-1]*(RC[14-1]*u[14-1]+RC[17-1]*u[1-1])+2.0*DJ[4-1]*u[17-1]+DJ[5-1]*u[16-1]
    s3 = s2+DJ[16-1]*(u[22-1]+u[28-1]+u[47-1]+u[49-1]+u[50-1]+u[52-1]+u[48-1]+u[53-1])
    s4 = s3+u[14-1]*(0.15*RC[123-1]*u[9-1]+0.8e-1*RC[160-1]*u[44-1])
    s5 = s4
    s7 = 0.55*RC[150-1]*u[14-1]*u[41-1]
    s9 = -1
    s12 = RC[222-1]*u[30-1]+RC[75-1]*u[12-1]+RC[81-1]*u[7-1]+RC[87-1]*u[52-1]+RC[86-1]* \
          u[13-1]+RC[235-1]*u[48-1]+RC[109-1]*u[8-1]+RC[125-1]*u[9-1]+RC[234-1]*u[10-1]+ \
          RC[223-1]*u[53-1]+RC[219-1]*u[34-1]+RC[31-1]*u[17-1]+RC[21-1]*u[2-1]+RC[148-1]* \
          u[62-1]+RC[64-1]*u[20-1]+RC[39-1]*u[3-1]+RC[71-1]*u[6-1]
    s11 = s12+RC[30-1]*u[15-1]+RC[59-1]*u[5-1]+RC[70-1]*u[4-1]+RC[35-1]*u[16-1]+RC[13-1]* \
          u[14-1]+RC[221-1]*u[32-1]+RC[68-1]*(u[22-1]+u[28-1]+u[47-1]+u[50-1]+u[51-1]+u[49-1])+ \
          RC[66-1]*u[11-1]+RC[151-1]*u[41-1]+RC[153-1]*u[44-1]+RC[63-1]*u[46-1]+RC[33-1]*u[18-1]+ \
          RC[158-1]*u[54-1]+RC[146-1]*u[63-1]+RC[149-1]*(u[65-1]+u[66-1])+RC[147-1]*u[64-1]+ \
          RC[161-1]*u[55-1]
    s12 = u[37-1]
    s10 = s11*s12
    s8 = s9*s10
    s6 = s7+s8
    
    deriv = [DJ[3-1]*u[2-1]+DJ[13-1]*u[39-1]+RC[19-1]*u[2-1]*u[39-1]+EMIS1/HMIX- \
             (RC[5-1]*u[36-1]+RC[11-1]*u[14-1]+RC[17-1]*u[15-1]+RC[72-1]*u[23-1]+RC[79-1]* \
              (u[24-1]+ u[57-1])+RC[15-1]*u[39-1]+RC[60-1]*(u[19-1]+u[26-1]+u[27-1]+u[29-1]+u[31-1]+u[33-1]+ \
                                                            u[35-1]+u[43-1]+u[45-1]+u[59-1]+u[61-1]+u[60-1]))*u[1-1], \

             u[1-1]*(RC[5-1]*u[36-1]+RC[11-1]*u[14-1]+RC[17-1]*u[15-1]+RC[72-1]*u[23-1]+ \
                     RC[79-1]*(u[24-1]+u[57-1])+0.2e+1*RC[15-1]*u[39-1])+RC[60-1]*u[1-1]* \
             (u[19-1]+u[26-1]+u[27-1]+u[29-1]+u[31-1]+u[33-1]+u[35-1]+u[59-1])+ \
             RC[60-1]*u[1-1]*(0.86*u[43-1]+0.19e+1*u[61-1]+0.11e+1*u[60-1]+0.95*u[45-1])+ \
             DJ[14-1]*u[39-1]+DJ[5-1]*u[16-1]+DJ[15-1]*u[40-1]+RC[29-1]*u[40-1]+RC[78-1]* \
             (u[25-1]+u[58-1])-(DJ[3-1]+RC[12-1]*u[14-1]+RC[20-1]*u[39-1]+RC[21-1]*u[37-1]+RC[48-1]+RC[77-1]* \
                                (u[24-1]+u[57-1]))*u[2-1], \

             EMIS3/HMIX-(RC[39-1]*u[37-1]+RC[40-1]*u[19-1]+RC[47-1])*u[3-1], \

             EMIS4/HMIX+u[37-1]*(RC[66-1]*u[11-1]+2.0*RC[221-1]*u[32-1]+RC[222-1]* u[30-1])+ \
             u[14-1]*(0.44*RC[112-1]*u[8-1]+0.4*RC[123-1]*u[9-1]+0.5e-1*RC[160-1]*u[44-1]+0.5e-1* \
                      RC[150-1]*u[41-1])+u[11-1]*(DJ[6-1]+DJ[7-1]+RC[69-1]*u[39-1])+DJ[8-1]*u[12-1]+DJ[11-1]* \
             u[30-1]+2.0*DJ[7-1]*u[32-1]-RC[70-1]*u[37-1]*u[4-1], \

             EMIS5/HMIX+0.7e-1*RC[123-1]*u[14-1]*u[9-1]-RC[59-1]*u[37-1]*u[5-1], \

             EMIS6/HMIX-RC[71-1]*u[37-1]*u[6-1], \

             EMIS7/HMIX-RC[81-1]*u[37-1]*u[7-1], \

             EMIS8/HMIX-(RC[109-1]*u[37-1]+RC[112-1]*u[14-1])*u[8-1], \

             EMIS9/HMIX+0.7e-1*RC[150-1]*u[14-1]*u[41-1]-(RC[123-1]*u[14-1]+RC[125-1]*u[37-1])*u[9-1], \

             EMIS10/HMIX-RC[234-1]*u[37-1]*u[10-1], \

             u[19-1]*(RC[60-1]*u[1-1]+(2.0*RC[61-1]+RC[62-1])*u[19-1]+RC[80-1]*u[24-1]+RC[40-1]*u[3-1])+ \
             u[37-1]*(RC[63-1]*u[46-1]+RC[67-1]*u[22-1])+u[1-1]*RC[60-1]* \
             (2.0*u[29-1]+u[31-1]+0.74*u[43-1]+0.266*u[45-1]+0.15*u[60-1])+u[14-1]* \
             (0.5*RC[123-1]*u[9-1]+RC[112-1]*u[8-1]+0.7*RC[157-1]*u[56-1]+0.8*RC[160-1]*u[44-1]+0.8*RC[150-1]* \
              u[41-1])+2.0*DJ[7-1]*u[32-1]+DJ[16-1]* \
             (u[22-1]+0.156e+1*u[50-1]+u[51-1])-(RC[66-1]*u[37-1]+DJ[6-1]+DJ[7-1]+RC[69-1]*u[39-1]+RC[53-1])*u[11-1], \

             u[1-1]*(RC[72-1]*u[23-1]+RC[83-1]*u[26-1]*RPATH4+RC[105-1]*u[27-1]+RC[126-1]*u[31-1]+ \
                     0.95*RC[162-1]*u[61-1]+0.684*RC[154-1]*u[45-1])+u[37-1]* \
             (RC[64-1]*u[20-1]+RC[76-1]*u[28-1]+RC[76-1]*u[50-1])+0.5*RC[123-1]*u[14-1]*u[9-1]+0.4e-1* \
             RC[160-1]*u[14-1]*u[44-1]+DJ[16-1]*(u[28-1]+0.22*u[50-1]+0.35*u[49-1]+u[51-1]+u[52-1])- \
             (DJ[8-1]+RC[75-1]*u[37-1]+RC[53-1])*u[12-1], \

             RC[83-1]*u[1-1]*u[26-1]*RPATH3+(0.65*DJ[16-1]+RC[76-1]*u[37-1])*u[49-1]+RC[76-1]*u[37-1]*u[51-1]+ \
             (RC[159-1]*u[1-1]+DJ[16-1])*u[59-1]+0.95*RC[162-1]*u[1-1]*u[61-1]-(DJ[9-1]+RC[86-1]*u[37-1]+RC[53-1])*u[13-1], \

             RC[1-1]*u[36-1]+RC[89-1]*u[15-1]*u[24-1]-(RC[11-1]*u[1-1]+RC[12-1]*u[2-1]+RC[13-1]*u[37-1]+RC[14-1]*u[15-1]+RC[49-1]+ \
                                                       RC[112-1]*u[8-1]+RC[123-1]*u[9-1]+RC[157-1]*u[56-1]+RC[160-1]*u[44-1]+RC[150-1]* \
                                                       u[41-1]+DJ[1-1]+DJ[2-1])*u[14-1], \

             s1+2.0*DJ[6-1]*u[11-1]+DJ[16-1]*(u[22-1]+u[28-1]+0.65*u[49-1]+u[50-1]+u[51-1]+u[48-1]+u[53-1])+ \
             u[39-1]*(RC[26-1]*u[17-1]+RC[69-1]*u[11-1])+u[14-1]*(0.12*RC[112-1]*u[8-1]+0.28*RC[123-1]*u[9-1]+ \
                                                                  0.6e-1*RC[160-1]*u[44-1])+0.6e-1*RC[150-1]*u[14-1]*u[41-1]- \
             (RC[14-1]*u[14-1]+RC[17-1]*u[1-1]+RC[30-1]*u[37-1]+2.0*RC[36-1]* \
              u[15-1]+RC[65-1]*u[19-1]+RC[74-1]*u[23-1]+(RC[88-1]+RC[89-1])*u[24-1]+ \
              RC[85-1]*(u[26-1]+u[29-1]+u[31-1]+u[27-1]+u[57-1]+u[45-1]+u[61-1]+u[59-1]+u[33-1]+ \
                        u[35-1]+u[43-1]+u[60-1]))*u[15-1], \

             RC[21-1]*u[2-1]*u[37-1]+u[39-1]*(RC[26-1]*u[17-1]+RC[69-1]*u[11-1])-(RC[35-1]*u[37-1]+DJ[5-1]+RC[45-1])*u[16-1], \

             RC[36-1]*u[15-1]**2-(RC[31-1]*u[37-1]+DJ[4-1]+RC[43-1]+RC[26-1]*u[39-1]+RC[47-1])*u[17-1], \

             DJ[7-1]*u[11-1]+u[14-1]*(0.13*RC[112-1]*u[8-1]+0.7e-1*RC[123-1]*u[9-1])-RC[33-1]*u[37-1]*u[18-1], \

             u[37-1]*(RC[59-1]*u[5-1]+RC[68-1]*u[22-1])+u[24-1]*(RC[79-1]*u[1-1]+2.0*RC[94-1]*u[24-1])+ \
             DJ[8-1]*u[12-1]+DJ[16-1]*u[47-1]+0.31*RC[123-1]*u[14-1]*u[9-1]-(RC[40-1]*u[3-1]+RC[60-1]*u[1-1]+ \
                                                                             2.0*RC[61-1]*u[19-1]+2.0*RC[62-1]*u[19-1]+ \
                                                                             RC[65-1]*u[15-1]+0.5*RC[80-1]*u[24-1])*u[19-1], \

             EMIS13/HMIX-RC[64-1]*u[37-1]*u[20-1], \

             (RC[40-1]*u[19-1]+RC[39-1]*u[37-1])*u[3-1]+0.5e-1*EMIS3/HMIX-RC[51-1]*u[21-1], \

             RC[65-1]*u[15-1]*u[19-1]-(RC[43-1]+DJ[16-1]+(RC[67-1]+RC[68-1])*u[37-1])*u[22-1], \

             u[37-1]*(RC[71-1]*u[6-1]+RC[68-1]*u[28-1])+0.35*DJ[16-1]*u[49-1]+RC[83-1]*u[1-1]*u[26-1]* \
             RPATH4+DJ[9-1]*u[13-1]-(RC[72-1]*u[1-1]+RC[74-1]*u[15-1])*u[23-1], \

             u[37-1]*(RC[75-1]*u[12-1]+RC[222-1]*u[30-1]+RC[68-1]*u[47-1])+RC[105-1]* \
             u[1-1]*u[27-1]+RC[78-1]*u[25-1]+DJ[11-1]*u[30-1]+DJ[9-1]*u[13-1]+DJ[16-1]*u[52-1]+ \
             0.684*RC[154-1]*u[1-1]*u[45-1]-(RC[77-1]*u[2-1]+RC[79-1]*u[1-1]+RC[80-1]*u[19-1]+ \
                                             2.0*RC[94-1]*u[24-1]+(RC[88-1]+RC[89-1])*u[15-1])*u[24-1], \

             RC[77-1]*u[24-1]*u[2-1]-(RC[50-1]+RC[78-1])*u[25-1], \

             u[37-1]*(RC[81-1]*u[7-1]+RC[68-1]*u[49-1])-(RC[83-1]*u[1-1]+RC[85-1]*u[15-1])*u[26-1], \

             u[37-1]*(RC[86-1]*u[13-1]+RC[87-1]*u[52-1])-(RC[105-1]*u[1-1]+RC[85-1]*u[15-1])*u[27-1], \

             RC[74-1]*u[15-1]*u[23-1]-((RC[76-1]+RC[68-1])*u[37-1]+DJ[16-1]+RC[52-1])*u[28-1], \

             u[37-1]*(RC[109-1]*u[8-1]+RC[68-1]*u[50-1])-(RC[110-1]*u[1-1]+RC[85-1]*u[15-1])*u[29-1], \

             RC[236-1]*u[1-1]*u[33-1]+RC[220-1]*u[1-1]*u[35-1]+0.266*RC[154-1]*u[1-1]*u[45-1]+ \
             0.82*RC[160-1]*u[14-1]*u[44-1]+DJ[16-1]*(u[48-1]+u[53-1])-(DJ[11-1]+RC[222-1]*u[37-1])*u[30-1], \

             u[37-1]*(RC[125-1]*u[9-1]+RC[68-1]*u[51-1])-(RC[126-1]*u[1-1]+RC[85-1]*u[15-1])*u[31-1], \

             RC[220-1]*u[1-1]*u[35-1]+DJ[16-1]*u[53-1]-(2.0*DJ[7-1]+RC[221-1]*u[37-1])*u[32-1], \

             u[37-1]*(RC[234-1]*u[10-1]+RC[235-1]*u[48-1])-(RC[236-1]*u[1-1]+RC[85-1]*u[15-1])*u[33-1], \

             RC[236-1]*u[1-1]*u[33-1]+DJ[16-1]*u[48-1]-RC[219-1]*u[37-1]*u[34-1], \

             u[37-1]*(RC[219-1]*u[34-1]+RC[223-1]*u[53-1])-(RC[220-1]*u[1-1]+RC[85-1]*u[15-1])*u[35-1], \

             DJ[1-1]*u[14-1]+DJ[3-1]*u[2-1]+DJ[14-1]*u[39-1]+RC[7-1]*u[38-1]+0.2*RC[160-1]*u[14-1]* \
             u[44-1]+0.3*RC[150-1]*u[14-1]*u[41-1]-(RC[1-1]+RC[5-1]*u[1-1])*u[36-1], \

             s5+s6, \

             DJ[2-1]*u[14-1]-(RC[7-1]+RC[8-1]*H2O)*u[38-1], \

             (RC[29-1]+DJ[15-1])*u[40-1]+RC[12-1]*u[14-1]*u[2-1]+RC[35-1]*u[37-1]*u[16-1]- \
             (RC[15-1]*u[1-1]+RC[26-1]*u[17-1]+RC[163-1]*u[41-1]+RC[19-1]*u[2-1]+RC[20-1]* \
              u[2-1]+DJ[13-1]+DJ[14-1]+RC[69-1]*u[11-1])*u[39-1], \

             RC[20-1]*u[39-1]*u[2-1]-(RC[29-1]+DJ[15-1]+RC[45-1])*u[40-1], \

             EMIS11/HMIX-(RC[151-1]*u[37-1]+RC[163-1]*u[39-1]+RC[150-1]*u[14-1])*u[41-1], \

             RC[45-1]*u[16-1]+2.0*RC[44-1]*u[40-1]-RC[51-1]*u[42-1], \

             u[37-1]*(RC[151-1]*u[41-1]+RC[156-1]*u[56-1])+0.12*RC[152-1]*u[1-1]* \
             u[43-1]-(RC[152-1]*u[1-1]+RC[155-1]*u[15-1])*u[43-1], \

             RC[60-1]*u[1-1]*(0.42*u[43-1]+0.5e-1*u[60-1])+0.26*RC[150-1]*u[14-1]* \
             u[41-1]-(RC[153-1]*u[37-1]+RC[160-1]*u[14-1])*u[44-1], \

             RC[153-1]*u[44-1]*u[37-1]+RC[148-1]*u[37-1]*u[62-1]-(RC[154-1]*u[1-1]+RC[85-1]*u[15-1])*u[45-1], \

             RC[62-1]*u[19-1]**2-RC[63-1]*u[37-1]*u[46-1], \

             RC[88-1]*u[15-1]*u[24-1]-(RC[68-1]*u[37-1]+DJ[16-1]+RC[52-1])*u[47-1], \

             RC[85-1]*u[15-1]*u[33-1]-(RC[235-1]*u[37-1]+DJ[16-1]+RC[52-1])*u[48-1], \

             RC[85-1]*u[15-1]*u[26-1]-((RC[76-1]+RC[68-1])*u[37-1]+DJ[16-1]+RC[52-1])*u[49-1], \

             RC[85-1]*u[15-1]*u[29-1]-((RC[76-1]+RC[68-1])*u[37-1]+DJ[16-1]+RC[52-1])*u[50-1], \

             RC[85-1]*u[15-1]*u[31-1]-((RC[76-1]+RC[68-1])*u[37-1]+DJ[16-1]+RC[52-1])*u[51-1], \

             RC[85-1]*u[15-1]*u[27-1]-(RC[87-1]*u[37-1]+DJ[16-1]+RC[52-1])*u[52-1], \

             RC[85-1]*u[15-1]*u[35-1]-(RC[223-1]*u[37-1]+DJ[16-1]+RC[52-1])*u[53-1], \

             RC[60-1]*u[1-1]*(0.32*u[43-1]+0.1*u[60-1])+0.67*RC[150-1]*u[14-1]*u[41-1]-RC[158-1]*u[37-1]*u[54-1], \

             RC[60-1]*u[1-1]*(0.14*u[43-1]+0.5e-1*u[45-1]+0.85*u[60-1])-RC[161-1]*u[37-1]*u[55-1], \

             RC[155-1]*u[15-1]*u[43-1]-(RC[156-1]*u[37-1]+RC[157-1]*u[14-1]+RC[52-1])*u[56-1], \

             0.5*RC[158-1]*u[37-1]*u[54-1]+RC[78-1]*u[58-1]+RC[149-1]*u[37-1]*u[66-1]- \
             (RC[77-1]*u[2-1]+RC[79-1]*u[1-1]+RC[85-1]*u[15-1])*u[57-1], \

             RC[77-1]*u[57-1]*u[2-1]-(RC[50-1]+RC[78-1])*u[58-1], \

             RC[79-1]*u[1-1]*u[57-1]+RC[146-1]*u[37-1]*u[63-1]-(RC[159-1]*u[1-1]+RC[85-1]*u[15-1])*u[59-1], \

             RC[163-1]*u[39-1]*u[41-1]+RC[147-1]*u[37-1]*u[64-1]-(RC[164-1]*u[1-1]+RC[85-1]*u[15-1])*u[60-1], \

             RC[161-1]*u[37-1]*u[55-1]+RC[149-1]*u[37-1]*u[65-1]-(RC[162-1]*u[1-1]+RC[85-1]*u[15-1])*u[61-1], \

             RC[85-1]*u[15-1]*u[45-1]-(RC[148-1]*u[37-1]+RC[52-1])*u[62-1], \

             RC[85-1]*u[15-1]*u[59-1]-(RC[146-1]*u[37-1]+RC[52-1])*u[63-1], \

             RC[85-1]*u[15-1]*u[60-1]-(RC[147-1]*u[37-1]+RC[52-1])*u[64-1], \

             RC[85-1]*u[15-1]*u[61-1]-(RC[149-1]*u[37-1]+RC[52-1])*u[65-1], \

             RC[85-1]*u[15-1]*u[57-1]-(RC[149-1]*u[37-1]+RC[52-1])*u[66-1]]
    return np.array(deriv)

t = tt[0]
dt = 1.0
tn = tt[-1] 
# initial conditions
u = np.empty(66)
for i in range(0,13):
    u[i] = 1.e+7
for i in range(13, 66):
    u[i] = 100.
u[0] = 1.e+9
u[1] = 5.e+9
u[2] = 5.e+9
u[3] = 3.8e+12
u[4] = 3.5e+13
u[13] = 5.e+11
u[37] = 1.e-3

