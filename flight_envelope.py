########################################################
#Script to plot flight envelopes and extract data 
########################################################

#imports
import argparse, numpy, pdb, sys, os
import matplotlib.pyplot as plt
sys.path.append('../utils')
from read_data import read_dat
from standard_atmosphere import density_imperial,density_SI

#Functions
def gust_load_factor(W,S,rho,c,a,g,V,Ude):
    '''Calculates the gust load factor according to FAR 23 Sec.23.341(c)'''
    mi_g = 2*(W/S)/(rho*c*a*g)
    Kg = 0.88*mi_g/(5.3 + mi_g)
    n = 1 + Kg*Ude*V*a/(489*(W/S))
    return n
    
def derived_gust_velocity(V_type,h):
    '''Calculates Ude according to FAR 23 Sec.23.333(c)(1), h in ft, Ud in ft/s'''
    if V_type == 'VC':
        if h < 20000:
            Ude = 50
        elif h < 50000:
            Ude = 50 - (50-25)*(h-20000)/(50000-20000)
    elif V_type == 'VD':
        if h < 20000:
            Ude = 25
        elif h < 50000:
            Ude = 25 - (25-12.5)*(h-20000)/(50000-20000)
    elif V_type == 'VB':
        if h < 20000:
            Ude = 66
        elif h < 50000:
            Ude = 66 - (66-38)*(h-20000)/(50000-20000)
    elif V_type == 'Vf':
        Ude = 25
    return Ude

def derived_gust_velocity_FAR25(V_type,h):
    '''Calculates Uds according to FAR 25, h in ft, Ud in ft/s'''
    if V_type == 'VC':
        if h < 15000:
            Ude = 56 - (56-44)*(h)/(15000)
        elif h < 50000:
            Ude = 44 - (44-26)*(h-15000)/(50000-15000)
    elif V_type == 'VD':
        if h < 15000:
            Ude = 28
            Ude = 28 - (28-22)*(h)/(15000)
        elif h < 50000:
            Ude = 28 - (28-22)*(h-15000)/(50000-15000)
    elif V_type == 'Vf':
        Ude = 25
    return Ude

def limit_maneuvering_load_factor(cat,W,n):
    '''Return the positive and negative limit load factors according to 
    FAR 23 Sec.23.337 and FAR 25. If the load factor specified in the input is less then
    the minimun required, it raises an error.'''
    if cat == 'FAR25':
        n_pos = numpy.ceil((2.1 + 24000/(W+10000))*10)/10
        if n_pos <= 2.5:
            n_pos = 2.5
        n_neg = 1.0
    if cat in ['normal','commuter']:
        n_pos = numpy.ceil((2.1 + 24000/(W+10000))*10)/10
        if n_pos > 3.8:
            n_pos = 3.8
        if n_pos > n:
            sys.exit('Specified load factor is less then the minimum required of {}'.format(n_pos))
        elif n > 3.8:
            print 'Specified load factor doesn\'t need to be higher than 3.8. Calculations will continue with n={}'.format(n)
        n_pos = n
        n_neg = numpy.ceil(4*n_pos)/10
    if cat == 'utility':
        n_pos = 4.4
        if n_pos > n:
            sys.exit('Specified load factor is less then the minimum required of {}'.format(n_pos))
        n_pos = n
        n_neg = numpy.ceil(4*n_pos)/10
    if cat == 'acrobatic':
        n_pos = 6.0
        if n_pos > n:
            sys.exit('Specified load factor is less then the minimum required of {}'.format(n_pos))
        n_pos = n
        n_neg = numpy.ceil(5*n_pos)/10
    return n_pos, n_neg

def extract_velocities(cat,cl_max,cl_max_f,cl_min,W,S,rho,n_pos,n_neg):
    '''Get Vs, VA, VC and VD according to FAR 23 Sec.23.335 and FAR 25 Sec.25.335'''
    #VB will be considered equal VA
    Vs_pos = numpy.ceil(numpy.sqrt(W/(cl_max*0.5*rho*S))*0.514444) #to convert to knots
    Vsf = numpy.ceil(numpy.sqrt(W/(cl_max_f*0.5*rho*S))*0.514444) #to convert to knots
    Vs_neg = numpy.ceil(numpy.sqrt(W/(numpy.absolute(cl_min)*0.5*rho*S))*0.514444) #to convert to knots
    VA_pos = numpy.ceil(Vs_pos*numpy.sqrt(n_pos)) 
    VA_neg = numpy.ceil(Vs_neg*numpy.sqrt(n_neg)) 
    if cat in ['normal','commuter']:
        if W/S < 20:
            VC = numpy.ceil(33*numpy.sqrt(W/S))
            VD = numpy.ceil(1.4*VC)
        else:
            VC = numpy.ceil((33-(33-28.6)*(W/S-20)/(100-20))*numpy.sqrt(W/S))
            VD = numpy.ceil((1.4-(1.4-1.35)*(W/S-20)/(100-20))*VC)
    elif cat == 'FAR25':
        VC = float(input('Design cruise speed(m/s):'))*0.514444
        VD = 1/0.8*VC
        print 'VD should be in the interval [{},{}] m/s.'.format(VC/0.514444,VD/0.514444)
        VD = float(input('Specify the design dive speed(m/s):'))*0.514444
    elif cat == 'utility':
        if W/S < 20:
            VC = numpy.ceil(33*numpy.sqrt(W/S))
            VD = numpy.ceil(1.5*VC)
        else:
            VC = numpy.ceil((33-(33-28.6)*(W/S-20)/(100-20))*numpy.sqrt(W/S))
            VD = numpy.ceil((1.5-(1.5-1.35)*(W/S-20)/(100-20))*VC)
    elif cat == 'acrobatic':
        if W/S < 20:
            VC = numpy.ceil(36*numpy.sqrt(W/S))
            VD = numpy.ceil(1.55*VC)
        else:
            VC = numpy.ceil((36-(36-28.6)*(W/S-20)/(100-20))*numpy.sqrt(W/S))
            VD = numpy.ceil((1.55-(1.55-1.35)*(W/S-20)/(100-20))*VC)
    if VC < Vs_pos or VC < Vs_neg:
        sys.exit('Stall speed higher than minimum cruise speed, please reconsider your design.')
    if VC < VA_pos:
        VC = numpy.ceil(1.1*VA_pos)
    if VD < 1.25*VC:
        VD = numpy.ceil(1.25*VC) 
    return Vs_pos, Vs_neg, Vsf, VA_pos, VA_neg, VC, VD

def write_log(data_dict,Vs_pos,Vs_neg,Vf,Vsf,VA_pos,VA_neg,VC,VD,ac_dir,plot_dict,tot_dict,n_pos,n_neg):
    '''Write the relevant data in a log file, based on the template.'''

    a = float(data_dict['a'][0])
    output = data_dict['output'][0]
    h = float(data_dict['h'][0])
    ac_name = ' '.join(x for x in data_dict['ac_name'])
    W = float(data_dict['W'][0])
    cat = data_dict['cat'][0]
    n = float(data_dict['n'][0])
    cl_max_f = float(data_dict['cl_max_f'][0])
    cl_max = float(data_dict['cl_max'][0])
    cl_min = float(data_dict['cl_min'][0])
    c = float(data_dict['c'][0])
    S = float(data_dict['S'][0])
    g = float(data_dict['g'][0])
    if output == 'imperial':
        rho = density_imperial(h)
        rho_un = 'sl/ft3'
        vel_un = ' knots'
    else:
        rho = density_SI(h*0.3048)
        rho_un = ' kg/m3'
        vel_un = '   m/s'
    W = str(W)
    while len(W)<7:
        W = ' '+W
    c = str(c)
    while len(c)<7:
        c = ' '+c
    S = str(S)
    while len(S)<7:
        S = ' '+S
    h = str(h)
    while len(h)<7:
        h = ' '+h
    a = str(a)
    while len(a)<7:
        a = ' '+a
    g = str(g)
    while len(g)<7:
        g = ' '+g
    n = str(n)
    while len(n)<4:
        n = ' '+n
    try:
        rho = str(rho)[0:6]
    except:
        rho = str(rho)
    while len(rho)<7:
        rho = ' '+rho
    cl_max = str(cl_max)
    while len(cl_max)<7:
        cl_max = ' '+cl_max
    cl_min = str(cl_min)
    while len(cl_min)<9:
        cl_min = ' '+cl_min
    cl_max_f = str(cl_max_f)
    while len(cl_max_f)<9:
        cl_max_f = ' '+cl_max_f
    while len(cl_max_f)<8:
        cl_max_f = ' '+cl_max_f
    while len(output)<8:
        output = ' '+output
    W_unit = data_dict['W'][1]
    while len(W_unit) < 6:
        W_unit = ' '+W_unit
    S_unit = data_dict['S'][1]
    while len(S_unit) < 6:
        S_unit = ' '+S_unit
    c_unit = data_dict['c'][1]
    while len(c_unit) < 6:
        c_unit = ' '+c_unit
    h_unit = data_dict['h'][1]
    while len(h_unit) < 6:
        h_unit = ' '+h_unit
    g_unit = data_dict['g'][1]
    while len(g_unit) < 6:
        g_unit = ' '+g_unit
    Vf_l_ac = str(tot_dict['pos'][1][int(Vf*10)])[0:4]
    while len(Vf_l_ac) < 7:
        Vf_l_ac = ' '+Vf_l_ac
    Vf_u_ac = str(1.5*tot_dict['pos'][1][int(Vf*10)])[0:4]
    while len(Vf_u_ac) < 7:
        Vf_u_ac = ' '+Vf_u_ac
    VA_pos_l_ac = str(tot_dict['pos'][1][int(VA_pos*10)])[0:4]
    while len(VA_pos_l_ac) < 11:
        VA_pos_l_ac = ' '+VA_pos_l_ac
    VA_pos_u_ac = str(1.5*tot_dict['pos'][1][int(VA_pos*10)])[0:4]
    while len(VA_pos_u_ac) < 11:
        VA_pos_u_ac = ' '+VA_pos_u_ac
    VA_neg_l_ac = str(tot_dict['neg'][1][int(VA_neg*10)])[0:5]
    while len(VA_neg_l_ac) < 11:
        VA_neg_l_ac = ' '+VA_neg_l_ac
    VA_neg_u_ac = str(1.5*tot_dict['neg'][1][int(VA_neg*10)])[0:5]
    while len(VA_neg_u_ac) < 11:
        VA_neg_u_ac = ' '+VA_neg_u_ac
    VC_pos_l_ac = str(tot_dict['pos'][1][int(VC*10)])[0:5]
    while len(VC_pos_l_ac) < 11:
        VC_pos_l_ac = ' '+VC_pos_l_ac
    VC_pos_u_ac = str(1.5*tot_dict['pos'][1][int(VC*10)])[0:5]
    while len(VC_pos_u_ac) < 11:
        VC_pos_u_ac = ' '+VC_pos_u_ac
    VC_neg_l_ac = str(tot_dict['neg'][1][int(VC*10)])[0:5]
    while len(VC_neg_l_ac) < 11:
        VC_neg_l_ac = ' '+VC_neg_l_ac
    VC_neg_u_ac = str(1.5*tot_dict['neg'][1][int(VC*10)])[0:5]
    while len(VC_neg_u_ac) < 11:
        VC_neg_u_ac = ' '+VC_neg_u_ac
    VD_neg_l_ac = str(tot_dict['neg'][1][int(VD*10)])[0:5]
    while len(VD_neg_l_ac) < 11:
        VD_neg_l_ac = ' '+VD_neg_l_ac
    VD_neg_u_ac = str(1.5*tot_dict['neg'][1][int(VD*10)])[0:5]
    while len(VD_neg_u_ac) < 11:
        VD_neg_u_ac = ' '+VD_neg_u_ac
    VD_pos_l_ac = str(tot_dict['pos'][1][int(VD*10)])[0:5]
    while len(VD_pos_l_ac) < 11:
        VD_pos_l_ac = ' '+VD_pos_l_ac
    VD_pos_u_ac = str(1.5*tot_dict['pos'][1][int(VD*10)])[0:5]
    while len(VD_pos_u_ac) < 11:
        VD_pos_u_ac = ' '+VD_pos_u_ac

    if output.strip() == 'imperial':
        Vs_pos = str(numpy.round(Vs_pos,2))
        while len(Vs_pos) < 9:
            Vs_pos = ' '+Vs_pos
        Vs_neg = str(numpy.round(Vs_neg,2))
        while len(Vs_neg) < 9:
            Vs_neg = ' '+Vs_neg
        Vf = str(numpy.round(Vf,2))
        while len(Vf) < 6:
            Vf = ' '+Vf
        Vsf = str(numpy.round(Vsf,2))
        while len(Vsf) < 6:
            Vsf= ' '+Vsf
        VA_pos = str(numpy.round(VA_pos,2))
        while len(VA_pos) < 9:
            VA_pos = ' '+VA_pos
        VA_neg = str(numpy.round(VA_neg,2))
        while len(VA_neg) < 9:
            VA_neg = ' '+VA_neg
        VC = str(numpy.round(VC,2))
        while len(VC) < 6:
            VC = ' '+VC
        VD = str(numpy.round(VD,2))
        while len(VD) < 6:
            VD = ' '+VD
    else:
        Vs_pos = str(numpy.round(Vs_pos/0.514444,2))
        while len(Vs_pos) < 9:
            Vs_pos = ' '+Vs_pos
        Vs_neg = str(numpy.round(Vs_neg/0.514444,2))
        while len(Vs_neg) < 9:
            Vs_neg = ' '+Vs_neg
        Vf = str(numpy.round(Vf/0.514444,2))
        while len(Vf) < 6:
            Vf = ' '+Vf
        Vsf = str(numpy.round(Vsf/0.514444,2))
        while len(Vsf) < 6:
            Vsf= ' '+Vsf
        VA_pos = str(numpy.round(VA_pos/0.514444,2))
        while len(VA_pos) < 9:
            VA_pos = ' '+VA_pos
        VA_neg = str(numpy.round(VA_neg/0.514444,2))
        while len(VA_neg) < 9:
            VA_neg = ' '+VA_neg
        VC = str(numpy.round(VC/0.514444,2))
        while len(VC) < 6:
            VC = ' '+VC
        VD = str(numpy.round(VD/0.514444,2))
        while len(VD) < 6:
            VD = ' '+VD
    while len(ac_name) < 16:
        ac_name += ' '
    while len(cat) < 12:
        cat += ' '
    n_pos = str(numpy.round(n_pos,2))
    while len(n_pos) < 5:
        n_pos = ' '+n_pos
    n_neg = str(numpy.round(-n_neg,2))
    while len(n_neg) < 5:
        n_neg = ' '+n_neg
    vf_pos_g = str(numpy.round(plot_dict['Vf_pos_gust'][1][-1],2))
    while len(vf_pos_g) < 8:
        vf_pos_g = ' '+vf_pos_g
    vc_pos_g = str(numpy.round(plot_dict['VC_pos_gust'][1][-1],2))
    while len(vc_pos_g) < 8:
        vc_pos_g = ' '+vc_pos_g
    vc_neg_g = str(numpy.round(plot_dict['VC_neg_gust'][1][-1],2))
    while len(vc_neg_g) < 8:
        vc_neg_g = ' '+vc_neg_g
    vd_pos_g = str(numpy.round(plot_dict['VD_pos_gust'][1][-1],2))
    while len(vd_pos_g) < 8:
        vd_pos_g = ' '+vd_pos_g
    vd_neg_g = str(numpy.round(plot_dict['VD_neg_gust'][1][-1],2))
    while len(vd_neg_g) < 8:
        vd_neg_g = ' '+vd_neg_g

    template_dict= {'W____ac':W,'c____ac':c,'S____ac':S,'h____ac':h,'rho__ac':rho,'g____ac':g,'a____ac':a,'cl_m_f_ac':cl_max_f,'cl_m_ac':cl_max,'cl_min_ac':cl_min,'n_ac':n,'W_unit':W_unit,'c_unit':c_unit,'S_unit':S_unit,'h_unit':h_unit,'rho_un':rho_un,'g_unit':g_unit,'Vs_pos_ac':Vs_pos,'Vs_neg_ac':Vs_neg,'Vsf_ac':Vsf,'Vf__ac':Vf,'VA_pos_ac':VA_pos,'VA_neg_ac':VA_neg,'VC__ac':VC,'VD__ac':VD,'ac_name_________':ac_name,'ac_cat______':cat,'n_pos':n_pos,'n_neg':n_neg,'vf_pos_g':vf_pos_g,'vc_pos_g':vc_pos_g,'vc_neg_g':vc_neg_g,'vd_pos_g':vd_pos_g,'vd_neg_g':vd_neg_g,'vel_un':vel_un,'__output':output,'Vf_l_ac':Vf_l_ac,'Vf_u_ac':Vf_u_ac,'VA_pos_l_ac':VA_pos_l_ac,'VA_pos_u_ac':VA_pos_u_ac,'VA_neg_l_ac':VA_neg_l_ac,'VA_neg_u_ac':VA_neg_u_ac,'VC_pos_l_ac':VC_pos_l_ac,'VC_pos_u_ac':VC_pos_u_ac,'VC_neg_l_ac':VC_neg_l_ac,'VC_neg_u_ac':VC_neg_u_ac,'VD_pos_l_ac':VD_pos_l_ac,'VD_pos_u_ac':VD_pos_u_ac,'VD_neg_l_ac':VD_neg_l_ac,'VD_neg_u_ac':VD_neg_u_ac}
    if cat.strip() == 'commuter':
        template_dict['VA_pos_commuter']='\nVA positive gust             {}\n'.format(str(plot_dict['VB_pos_gust'][1][-1])[0:3])
        template_dict['VA_neg_commuter']='VA positive gust            {}'.format(str(plot_dict['VB_neg_gust'][1][-1])[0:4])
    else:
        template_dict['VA_pos_commuter']=''
        template_dict['VA_neg_commuter']=''
    with open(ac_dir+'/'+ac_name.strip()+'.log','w') as log:
        with open('log_template.txt','r') as templ:
            for line in templ:
                for key in template_dict:
                    if key in line:
                        line = line.replace(key,template_dict[key])
                log.write(line)
    return
        
#Main
def plot_envelope(dat_path, outpath):
    '''Plots the envelope calling the other functions.'''
    data_dict = read_dat(dat_path)
    # extracting the data from the input file.
    if data_dict['h'][1] == 'ft':
        h = float(data_dict['h'][0])
    elif data_dict['h'][1] == 'm':
        h = float(data_dict['h'][0])/0.3048 #converting to ft
    else:
        sys.exit('Altitude unit {} not supported, convert it to m or ft.'.format(data_dict['h'][1]))
    ac_name = ' '.join(x for x in data_dict['ac_name'])
    if data_dict['W'][1] == 'lbs':
        W = float(data_dict['W'][0])
    elif data_dict['W'][1] == 'N':
        W = float(data_dict['W'][0])*0.224809 #converting to lbs
    else:
        sys.exit('Weight unit {} not supported, convert it to N or lbs.'.format(data_dict['W'][1]))
    cat = data_dict['cat'][0]
    n = float(data_dict['n'][0])
    cl_max_f = float(data_dict['cl_max_f'][0])
    cl_max = float(data_dict['cl_max'][0])
    cl_min = float(data_dict['cl_min'][0])
    if data_dict['c'][1] == 'ft':
        c = float(data_dict['c'][0])
    elif data_dict['c'][1] == 'm':
        c = float(data_dict['c'][0])/0.3048 #converting to ft
    else:
        sys.exit('Mean geometric chord unit {} not supported, convert it to m or ft.'.format(data_dict['c'][1]))
    if data_dict['S'][1] == 'ft^2':
        S = float(data_dict['S'][0])
    elif data_dict['S'][1] == 'm^2':
        S = float(data_dict['S'][0])/0.3048**2 #converting to ft^2
    else:
        sys.exit('Mean geometric chord unit {} not supported, convert it to m^2 or ft^2.'.format(data_dict['S'][1]))
    if data_dict['g'][1] == 'ft/s^2':
        g = float(data_dict['g'][0])
    elif data_dict['g'][1] == 'm/s^2':
        g = float(data_dict['g'][0])/0.3048 #converting to ft/s^2
    else:
        sys.exit('Gravity acceleration unit {} not supported, convert it to m/s^2 or ft/s^2.'.format(data_dict['g'][1]))
    a = float(data_dict['a'][0])
    output = data_dict['output'][0]
    # creating additional data:
    rho = density_imperial(h)
    rho_SI = density_SI(h*0.3048)
    
    n_pos, n_neg = limit_maneuvering_load_factor(cat,W,n)
        
    Vs_pos, Vs_neg, Vsf, VA_pos, VA_neg, VC, VD = extract_velocities(cat,cl_max,cl_max_f,cl_min,W/0.224809,S*0.3048**2,rho_SI,n_pos,n_neg)

    if cat == 'FAR25':
        Vf = numpy.ceil(numpy.max([1.6*Vs_pos,1.8*Vsf]))    
    else:
        Vf = numpy.ceil(numpy.max([1.4*Vs_pos,1.8*Vsf]))    

    # creating gust cases
    gust_dict = {}
    if cat == 'commuter':
        VB_pos = VA_pos
        VB_neg = VA_neg
        vel_dict = {'Vf_pos_gust':Vf, 'VB_pos_gust':VB_pos, 'VB_neg_gust':VB_neg, 'VC_pos_gust':VC, 'VC_neg_gust':VC ,'VD_pos_gust':VD, 'VD_neg_gust':VD}
        for vel in vel_dict:
            Ude = derived_gust_velocity(vel[0:2],h)
            if 'neg' in vel:
                n_gust = 1 - gust_load_factor(W,S,rho,c,a,g,vel_dict[vel],Ude)
            else:
                n_gust = gust_load_factor(W,S,rho,c,a,g,vel_dict[vel],Ude)
            gust_dict[vel] = [vel_dict[vel],n_gust]
    elif cat == 'FAR25':
        vel_dict = {'Vf_pos_gust':Vf, 'VC_pos_gust':VC, 'VC_neg_gust':VC ,'VD_pos_gust':VD, 'VD_neg_gust':VD}
        for vel in vel_dict:
            Ude = derived_gust_velocity_FAR25(vel[0:2],h)
            if 'neg' in vel:
                n_gust = 1 - gust_load_factor(W,S,rho,c,a,g,vel_dict[vel],Ude)
            else:
                n_gust = gust_load_factor(W,S,rho,c,a,g,vel_dict[vel],Ude)
            gust_dict[vel] = [vel_dict[vel],n_gust]
    else:
        vel_dict = {'Vf_pos_gust':Vf, 'VC_pos_gust':VC, 'VC_neg_gust':VC ,'VD_pos_gust':VD, 'VD_neg_gust':VD}
        for vel in vel_dict:
            Ude = derived_gust_velocity(vel[0:2],h)
            if 'neg' in vel:
                n_gust = 1 - gust_load_factor(W,S,rho,c,a,g,vel_dict[vel],Ude)
            else:
                n_gust = gust_load_factor(W,S,rho,c,a,g,vel_dict[vel],Ude)
            gust_dict[vel] = [vel_dict[vel],n_gust]

    # generating data to plot
    plot_dict = {}
    # gust plots
    for key in gust_dict:
        if 'pos' in key and 'gust' in key:
            plot_dict[key] = [numpy.linspace(0,gust_dict[key][0],vel_dict[key]*10+1),numpy.linspace(1,gust_dict[key][1],vel_dict[key]*10+1)]
        else:
            plot_dict[key] = [numpy.linspace(0,gust_dict[key][0],vel_dict[key]*10+1),numpy.linspace(0,gust_dict[key][1],vel_dict[key]*10+1)]
    # maneuvre plots

    x_range = numpy.linspace(0,VD,VD*10+1)
    n_man = numpy.zeros_like(x_range)
    A = numpy.reshape([VA_pos**2 , VA_pos, Vs_pos**2, Vs_pos],(2,2))
    [a,b] = numpy.linalg.solve(A,[n_pos,1])
    for i,x in enumerate(x_range):
        y = a*x**2 + b*x
        if y < n_pos:
            n_man[i] = y
        else:
            n_man[i] = n_pos
    plot_dict['pos_maneuvre'] = [x_range,n_man] 

    x_range = numpy.linspace(0,VD,VD*10+1)
    n_man = numpy.zeros_like(x_range)
    A = numpy.reshape([VA_neg**2 , VA_neg, Vs_neg**2, Vs_neg],(2,2))
    try:
        [a,b] = numpy.linalg.solve(A,[n_neg,1])
    except:
        b=0
        a=1/Vs_neg**2
    for i,x in enumerate(x_range):
        y = a*x**2 + b*x
        if y < n_neg:
            n_man[i] = -y
        else:
            n_man[i] = -n_neg
    plot_dict['neg_maneuvre'] = [x_range,n_man] 
    plot_dict['VD_maneuvre'] = [[VD,VD],[n_pos,-n_neg]]

    # flap curve
    x_range = numpy.linspace(0,Vf,Vf*10+1)
    n_fl = numpy.zeros_like(x_range)
    A = numpy.reshape([VA_pos**2 , VA_pos, Vsf**2, Vsf],(2,2))
    [a,b] = numpy.linalg.solve(A,[n_pos,1])
    for i,x in enumerate(x_range):
        y = a*x**2 + b*x
        if y < 2.0:
            n_fl[i] = y
        else:
            n_fl[i] = 2.0
    plot_dict['flap'] = [x_range,n_fl]
    
    # total envelopei
    ##Positive
    tot_dict = {}
    x_pos = []
    n_tot = []
    for i,x in enumerate(plot_dict['flap'][0]):
        if plot_dict['flap'][1][i] >= plot_dict['pos_maneuvre'][1][i] and plot_dict['flap'][1][i] <= 1:
            x_pos.append(x)
            n_tot.append(plot_dict['flap'][1][i])
        elif plot_dict['Vf_pos_gust'][1][-1] > plot_dict['flap'][1][-1] and plot_dict['Vf_pos_gust'][1][-1] > plot_dict['pos_maneuvre'][1][len(plot_dict['flap'][0])]:
            x_pos.append(plot_dict['Vf_pos_gust'][0][-1])
            n_tot.append(plot_dict['Vf_pos_gust'][1][-1])
        elif plot_dict['flap'][1][i] > plot_dict['pos_maneuvre'][1][i]:
            x_pos.append(x)
            n_tot.append(plot_dict['flap'][1][i])
        elif plot_dict['flap'][1][i] < plot_dict['pos_maneuvre'][1][i]:
            x_pos.append(x)
            n_tot.append(plot_dict['pos_maneuvre'][1][i])
    for x in numpy.arange(Vf,VA_pos,0.1):
        x_pos.append(x)
        n_tot.append(plot_dict['pos_maneuvre'][1][10*x])
    link_gust_pos_n = []
    for i in numpy.arange(VA_pos,VC,0.1):
        if cat == 'commuter':
            n_gp = plot_dict['VB_pos_gust'][1][-1]-(plot_dict['VB_pos_gust'][1][-1]-plot_dict['VC_pos_gust'][1][-1])*(i-VB_pos)/(VC-VB_pos)
            link_gust_pos_n.append(n_gp)
        else:
            link_gust_pos_n.append(plot_dict['VC_pos_gust'][1][10*i])
    for i in numpy.arange(VC,VD+0.1,0.1):
        n_gp = plot_dict['VC_pos_gust'][1][-1]-(plot_dict['VC_pos_gust'][1][-1]-plot_dict['VD_pos_gust'][1][-1])*(i-VC)/(VD-VC)
        link_gust_pos_n.append(n_gp)

    for i,x in enumerate(numpy.arange(VA_pos,VD+0.1,0.1)):
        x_pos.append(x)
        n_tot.append(numpy.max([link_gust_pos_n[i],n_pos]))
    x_pos.append(VD)
    n_tot.append(0)
    tot_dict['pos'] = [x_pos,n_tot]

    ##negative
    x_neg = []
    n_tot = []
    for x in numpy.arange(0,VA_neg,0.1):
        x_neg.append(x)
        n_tot.append(plot_dict['neg_maneuvre'][1][10*x])
    link_gust_neg_n = []
    for i in numpy.arange(VA_neg,VC,0.1):
        if cat == 'commuter':
            n_gn = plot_dict['VB_neg_gust'][1][-1]-(plot_dict['VB_neg_gust'][1][-1]-plot_dict['VC_neg_gust'][1][-1])*(i-VB_neg)/(VC-VB_neg)
            link_gust_neg_n.append(n_gn)
        else:
            link_gust_neg_n.append(plot_dict['VC_neg_gust'][1][10*i])
    for i in numpy.arange(VC,VD+0.1,0.1):
        n_gn = plot_dict['VC_neg_gust'][1][-1]-(plot_dict['VC_neg_gust'][1][-1]-plot_dict['VD_neg_gust'][1][-1])*(i-VC)/(VD-VC)
        link_gust_neg_n.append(n_gn)
    for i,x in enumerate(numpy.arange(VA_neg,VD+0.1,0.1)):
        x_neg.append(x)
        n_tot.append(numpy.min([link_gust_neg_n[i],-n_neg]))
    x_neg.append(VD)
    n_tot.append(0)
    tot_dict['neg'] = [x_neg,n_tot]

    ac_dir = os.path.join(outpath,ac_name)
    plot_dir = os.path.join(outpath,ac_name,'envelopes')
    if not os.path.isdir(ac_dir):
        os.mkdir(ac_dir)
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)
    # plot limit envelope
    plt.figure()
    plt.ylabel('Load factor')
    plt.title('{} gust and maneuvre flight envelope'.format(ac_name))
    if output == 'imperial':
        plt.xlabel('Velocity (knots)')
        for key in plot_dict:
            if 'gust' in key:
                plt.plot(plot_dict[key][0],plot_dict[key][1],'k.', markersize=1.5)
            else:
                plt.plot(plot_dict[key][0],plot_dict[key][1],'b')
            plt.xlim([0,1.1*VD])
            plt.ylim([-1.5*n_neg,1.5*n_pos])
    elif output == 'SI':
        plt.xlabel('Velocity (m/s)')
        for key in plot_dict:
            if 'gust' in key:
                plt.plot(numpy.array(plot_dict[key][0])/0.514444,plot_dict[key][1],'k.', markersize=1.5)
            else:
                plt.plot(numpy.array(plot_dict[key][0])/0.514444,plot_dict[key][1],'b')
            plt.xlim([0,1.1*VD/0.514444])
            plt.ylim([-1.5*n_neg,1.5*n_pos])
    plt.grid()
    plt.savefig(plot_dir+'/gust_maneuvre_limit_envelope.png')
    
    plt.figure()
    plt.ylabel('Load factor')
    plt.xlim([0,1.1*VD])
    plt.ylim([-1.5*n_neg,1.5*n_pos])
    plt.title('{} resultant flight envelope (limit condition)'.format(ac_name))
    if output == 'imperial':
        for key in tot_dict:
            plt.xlabel('Velocity (m/s)')
            plt.xlim([0,1.1*VD])
            plt.ylim([-1.5*n_neg,1.5*n_pos])
            plt.plot(tot_dict[key][0],tot_dict[key][1],'b',linewidth=3.0)
    elif output == 'SI':
        for key in tot_dict:
            plt.xlabel('Velocity (m/s)')
            plt.plot(numpy.array(tot_dict[key][0])/0.514444,tot_dict[key][1],'b',linewidth=3.0)
            plt.xlim([0,1.1*VD/0.514444])
            plt.ylim([-1.5*n_neg,1.5*n_pos])
    plt.grid()
    plt.savefig(plot_dir+'/resultant_limit_envelope.png')

    plt.figure()
    plt.ylabel('Load factor')
    plt.title('{} resultant flight envelope (ultimate condition)'.format(ac_name))
    if output == 'imperial':
        for key in tot_dict:
            plt.xlabel('Velocity (m/s)')
            plt.plot(tot_dict[key][0],1.5*numpy.array(tot_dict[key][1]),'b',linewidth=3.0)
    elif output == 'SI':
        for key in tot_dict:
            plt.xlabel('Velocity (m/s)')
            plt.plot(numpy.array(tot_dict[key][0])/0.514444,1.5*numpy.array(tot_dict[key][1]),'b',linewidth=3.0)
    plt.grid()
    plt.savefig(plot_dir+'/resultant_ultimate_envelope.png')

    write_log(data_dict,Vs_pos,Vs_neg,Vf,Vsf,VA_pos,VA_neg,VC,VD,ac_dir,plot_dict,tot_dict,n_pos,n_neg)

########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to plot flight envelopes and extract data.')
    parser.add_argument('-dat', nargs='+', default='/home/jh/Documents/scripts/projects/flight_envelope/data_flight_envelope.dat', help='Path to the dat file containing the input variables.')
    parser.add_argument('-o', nargs='+', default='.', help='Desired path to save the results.')
    args = parser.parse_args()

    plot_envelope(args.dat, args.o)
