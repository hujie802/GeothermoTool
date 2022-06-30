import PySimpleGUI as sg
from random import randint
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate,stats
# import operator

"""
    Another simple table created from Input Text Elements.  This demo adds the ability to "navigate" around the drawing using
    the arrow keys. The tab key works automatically, but the arrow keys are done in the code below.
"""
sg.theme('Lightblue 2')  # No excuse for gray windows
# Show a "splash" type message so the user doesn't give up waiting
sg.popup_quick_message('Hang on for a moment, this will take a bit to create....', auto_close=True, non_blocking=True)

#========================== Functions======================================
def corr_TC_HF(pdt,T0,unmixing_method,mixing_method,fluid):
    #create tempoary
    # print(pdt)
    # K_test = pdt['TCt\n(W/mK)']
    K_test = pdt.iloc[:,0]
    print(K_test)
    # phi_test = pdt['Phi\n(%)']/100
    phi_test = pdt.iloc[:, 1]/100
    # T_in_situ = pdt['T\n(C)']
    T_in_situ = pdt.iloc[:,2]
    # P_test = pdt['P\n(MPa)']
    P_test = pdt.iloc[:,3]
    # G_test = pdt['G\n(C/km)']
    G_test = pdt.iloc[:, 4]
    T0 = float(T0)

    TC_air = 26e-3
    TC_oil = 0.14 # Physical properties of rocks, 334
    TC_water = 0.66

    #calc TC for solid  (Unmixing)
    if unmixing_method == 'Square-root':
        K_solid = ((K_test**0.5 - phi_test*TC_air**0.5)/(1 - phi_test))**2 #TC for solid
    elif unmixing_method == 'Geometric':
        K_solid = (K_test/TC_air**phi_test)**(1/(1 - phi_test))  # TC for solid GEOMETRIC
    elif unmixing_method == 'Harmonic':
        K_solid = (1-phi_test)/(1/K_test-phi_test/TC_air)
    elif unmixing_method == 'Arithmetic':
        K_solid = (K_test-phi_test*TC_air)/(1-phi_test)
    else: # No Unmixing
        K_solid = K_test

    #Temperature Correction
    K_solid_corrected = ((T0 + 273)*1473/(1473-T0-273))*(K_solid-1.05)*(1/(T_in_situ+273) - 1/1473) + 1.05
    #determine the fluid
    if fluid == 'Water':
        TC_fluid = TC_water
    elif fluid == 'Gas':
        TC_fluid = TC_air
    else:
        TC_fluid = TC_oil
    #Mixing again
    if mixing_method == 'Square-root':
        K_c = (K_solid_corrected**0.5*(1-phi_test)+TC_fluid**0.5*phi_test)**2
    elif mixing_method == 'Geometric':
        K_c = K_solid_corrected**(1-phi_test)*TC_fluid**phi_test
    elif mixing_method == 'Harmonic':
        K_c = 1/(phi_test/TC_fluid+(1-phi_test)/K_solid_corrected)
    elif mixing_method == 'Arithmetic':
        K_c = phi_test*TC_fluid+(1-phi_test)*K_solid_corrected

    K_cp = K_c+0.0005*P_test

    #Calc heat flow
    HF_c = K_cp*G_test
    return K_solid,K_cp,HF_c

#+++++++++++++++Geothermal Energy Estinmating+++++++++++++++++
def gen_rand_values(min_v, mean_v, max_v, method, n):
    """
    mean_v = mean value, min_v = min value, max_v = max value
    n = resample numbers
    method: 0 for normal, 1 for triangular, 2 for unifrom, 3 for Constant
    """
    if method == 'Normal':
        sd = max_v / 2 - min_v / 2
        rand = np.random.normal(mean_v, sd, n)
        rand[rand < 0] = 0.001  # delete 0
        return rand
    elif method == 'Laplace':
        sd = max_v / 2 - min_v / 2
        rand = np.random.laplace(mean_v, sd, n)
        rand[rand < 0] = 0.001  # delete 0
        # print(rand)
        return rand
    elif method == 'Triangular':
        rand = np.random.triangular(min_v, mean_v, max_v, n)
        rand[rand < 0] = 0.001  # delete 0
        return rand
    elif method == 'Uniform':
        rand = np.random.uniform(min_v, max_v, n)
        rand[rand < 0] = 0.001  # delete 0
        return rand
    else:
        return np.linspace(mean_v, mean_v, n)


def thermal_energy(rhow, A, Ts, rhor, Cr, Cw, Tr, phi, h):
    Qr = A * h * (rhor * Cr * (1 - phi/100) * (Tr - Ts)) * 1e6  # change km to m
    Qw = A * h * (rhow * Cw * phi/100 * (Tr - Ts)) * 1e6  # change km to m
    Qt = Qr + Qw
    return Qt, Qr, Qw

import scipy.stats as stats
def plotter(ax, name, rand_values, min_v, mean_v, max_v, method):
    """
    ax = ax project
    name = the title of the figure
    rand_values = resampled values
    mean_v = mean value, min_v = min value, max_v = max value
    n = resample numbers
    method:  normal, triangular, unifrom,  for Constant
    """
    bins = 20
    # print(method)

    if method == 'Normal':
        sd = max_v / 2 - min_v / 2
        x = np.linspace(sd * (-4), sd * (4)) + mean_v
        y = stats.norm.pdf(x, mean_v, sd)
        y_max = y.max()
        ax.set_ylim([0, y_max * 1.15])
        ax.axhline(y_max * 1.05, color='k')
        rand_values_sort = np.sort(rand_values)
        itz = len(rand_values)
        rand_values_first_2_5 = rand_values_sort[int(itz * 0.025)]  # confident interval
        rand_values_last_2_5 = rand_values_sort[int(itz * 0.975)]
        ax.text(rand_values_first_2_5 / 2 + rand_values_last_2_5 / 2, y_max * 1.06, '95%', ha='center')
        ax.text(rand_values_first_2_5, y_max * 1.06, '2.5% ', ha='right')
        ax.text(rand_values_last_2_5, y_max * 1.06, ' 2.5%', ha='left')
        ax.plot(x, y)
        ax.hist(rand_values, density=True, bins=bins, color='tomato')
        ax.axvspan(rand_values_first_2_5, rand_values_last_2_5, alpha=0.3, color='red')
        ax.set_title(name)
    elif method == 'Laplace':
        sd = max_v / 2 - min_v / 2
        x = np.linspace(sd * (-4), sd * (4)) + mean_v
        y = stats.laplace.pdf(x, mean_v, sd)
        y_max = y.max()
        ax.set_ylim([0, y_max * 1.15])
        ax.axhline(y_max * 1.05, color='k')
        rand_values_sort = np.sort(rand_values)
        itz = len(rand_values)
        rand_values_first_2_5 = rand_values_sort[int(itz * 0.025)]  # confident interval
        rand_values_last_2_5 = rand_values_sort[int(itz * 0.975)]
        ax.text(rand_values_first_2_5 / 2 + rand_values_last_2_5 / 2, y_max * 1.06, '95%', ha='center')
        ax.text(rand_values_first_2_5, y_max * 1.06, '2.5% ', ha='right')
        ax.text(rand_values_last_2_5, y_max * 1.06, ' 2.5%', ha='left')
        ax.plot(x, y)
        ax.hist(rand_values, density=True, bins=bins, color='tomato')
        ax.axvspan(rand_values_first_2_5, rand_values_last_2_5, alpha=0.3, color='red')
        ax.set_title(name)
    elif method == 'Triangular':
        x = np.linspace(min_v, max_v, 100)
        scale = max_v - min_v
        c = (mean_v - min_v) / scale
        y = stats.triang.pdf(x, c, min_v, scale)
        y_max = y.max()
        ax.set_ylim([0, y_max * 1.15])
        ax.axhline(y_max * 1.05, color='k')
        rand_values_sort = np.sort(rand_values)
        itz = len(rand_values)
        rand_values_first_2_5 = rand_values_sort[int(itz * 0.025)]
        rand_values_last_2_5 = rand_values_sort[int(itz * 0.975)]
        ax.text(rand_values_first_2_5 / 2 + rand_values_last_2_5 / 2, y_max * 1.06, '95%', ha='center')
        ax.text(rand_values_first_2_5, y_max * 1.06, '2.5% ', ha='right')
        ax.text(rand_values_last_2_5, y_max * 1.06, ' 2.5%', ha='left')
        ax.plot(x, y)
        ax.hist(rand_values, density=True, bins=bins, color='tomato')
        ax.set_title(name)
        ax.axvspan(rand_values_first_2_5, rand_values_last_2_5, alpha=0.3, color='red')
    elif method == 'Uniform':
        x = np.linspace(min_v, max_v, 100)
        p = 1 / (max_v - min_v)
        y = np.linspace(p, p, 100)
        y_max = y.max()
        ax.set_ylim([0, y_max * 1.15])
        ax.axhline(y_max * 1.05, color='k')
        rand_values_sort = np.sort(rand_values)
        itz = len(rand_values)
        rand_values_first_2_5 = rand_values_sort[int(itz * 0.025)]
        rand_values_last_2_5 = rand_values_sort[int(itz * 0.975)]
        ax.text(rand_values_first_2_5 / 2 + rand_values_last_2_5 / 2, y_max * 1.06, '95%', ha='center')
        ax.text(rand_values_first_2_5, y_max * 1.06, '2.5% ', ha='right')
        ax.text(rand_values_last_2_5, y_max * 1.06, ' 2.5%', ha='left')
        ax.plot(x, y)
        ax.hist(rand_values, density=True, bins=bins, color='tomato')
        ax.set_title(name)
        ax.axvspan(rand_values_first_2_5, rand_values_last_2_5, alpha=0.3, color='red')
    else:
        x = np.linspace(mean_v, mean_v, 100)
        # print(x)
        y = np.linspace(0, 1, 100)
        ax.plot(x, y)
        # ax.hist(rand_values,density = True,bins = bins,color = 'tomato')
        ax.set_title(name)
#+++++++++++++++Geothermal Energy Estinmating+++++++++++++++++

#========================== Functions end======================================

#========================== Plot   ======================================

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib
matplotlib.use('Agg')
def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
    if canvas.children:
        for child in canvas.winfo_children():
            child.destroy()
    if canvas_toolbar.children:
        for child in canvas_toolbar.winfo_children():
            child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
    figure_canvas_agg.draw()
    toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
    toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)

class Toolbar(NavigationToolbar2Tk):
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)

#========================== Plot  end======================================

#===================================== GUI=================================
#++++++++++++++++++++++++++++++++++++TC tab++++++++++++++++++++++++++++++++
#TCt = tested TC, TCs = tested TC for solid, Phi = porosity, T = in-situ Temperature, TC_c = corrected TC, G= Gradiend, HF = heat flow
MAX_ROWS, MAX_COLS, COL_HEADINGS = 15, 8, ('TCt\n(W/(m·K))','Phi\n(%)','T\n(°C)','P\n(MPa)','G\n(°C/km)', 'TCs\n(W/(m·K))','TCc\n(W/(m·K))','HF\n(mW/m²')
initial_TC_data = [[2.5,15,100,0,30,0,0,0]]*MAX_ROWS
# header_list = COL_HEADINGS
headings = list(COL_HEADINGS)
#input simulated table
columns0 = [[sg.T(r, size=(2, 2))] + [sg.Input(initial_TC_data[r][c], justification='c',size= (10,2),
                                               key=(r, c)) for c in range(MAX_COLS)] for r in range(MAX_ROWS)]
initial_folder = None
if initial_folder is None:
    initial_folder = os.getcwd()
columns1 = [
        [sg.Button('Import TC')],
        [sg.Text('Unmixing Law', background_color=sg.DEFAULT_BACKGROUND_COLOR,
              justification='center', size=(20, 1))],
        [sg.Combo(values=('Square-root', 'Geometric', 'Harmonic','Arithmetic','None'),
                 default_value='Square-root', key='Unmixing_Law_Combo')],
        [sg.Text('Mixing Law', background_color=sg.DEFAULT_BACKGROUND_COLOR,
                 justification='center', size=(20, 1))],
        [sg.Combo(values=('Square-root', 'Geometric', 'Harmonic','Arithmetic'),
                 default_value='Square-root', key='Mixing_Law_Combo')],
        [sg.Text('Room temperature', background_color=sg.DEFAULT_BACKGROUND_COLOR,
              justification='l', size=(20, 1))],
        [sg.Input(20, key= 'T0',justification='r',size= (20,1))],
        [sg.Text('Fluid', background_color=sg.DEFAULT_BACKGROUND_COLOR,
              justification='center', size=(20, 1))],
        [sg.Combo(values=('Water', 'Gas', 'Oil'),
                 default_value='Water', key='Fluid_Combo', size=(12, 1))],
        [sg.Button('Calc TC')],
        [sg.FileSaveAs('Export TC',key='Export_TC', initial_folder=initial_folder,enable_events = True,disabled=True,file_types = (('MS EXCEL', '.xlsx'),))]  # TODO: better names
            ]
#TC correction
tab_TC_layout = \
        [
        [sg.Text('Correct Thermal Conductivity', font='Default 16')],
        [sg.Text(' ' * 4,key = 'header_TC',visible = True)] + [sg.Text(s, key=s,visible = True ,enable_events=True, font='Courier 12', size=(7, 2)) for i, s in enumerate(COL_HEADINGS)],
        [sg.Col(columns0,key = 'Col0',visible = True),
        sg.Table(values=initial_TC_data, headings=headings, def_col_width = 8,
                                # background_color='black',
                                auto_size_columns=False,
                                display_row_numbers=True,
                                justification='center',
                                num_rows=15,
                                # alternating_row_color='black',
                                key='Table_TC',
                                row_height=30,visible = False,size = [50,200]),
         sg.Col(columns1)],
        ]
#++++++++++++++++++++++++++++++++Logs tab++++++++++++++++++++++++++++++++
#well log analysis
heading_log = ['Depth (m)','Temperature (°C)']
initial_log = [[0,0]]
columns_logs = [[sg.Button('Import Logs')],
                [sg.Text('Number of step')],
                [sg.Slider(orientation='h', range=(10, 100),key='number_of_stp_logs',default_value=50)],
                [sg.Button('Plot Logs',disabled= True)],
                [sg.Button('Fit Logs',disabled= True)],
                [sg.Text("Segment number"),sg.Combo(values = [1,2,3,4],default_value=1, key='Number_of_segments')],
                [sg.Button('Fit all',disabled= True)],
                [sg.Text("Fitting information")],
                [sg.Multiline(size=(20,15), font='Courier 8',key='Fit_information',disabled = True)], # fit output information, disabled means no input
                # [sg.Output(size=(15,10))],
                # [sg.Button('Export Logs',disabled= True)]
                ]
figure_w, figure_h = 400, 300
col_canvas = sg.Col([[sg.Canvas(size=(figure_w, figure_h), key='canvas_logs')],[sg.Canvas(key='controls_cv')]]) # cavas for logs
tab_Well_log_layout =[[sg.T('Temperature Logs',font='Default 18')],
                      [sg.Table(values=initial_log, headings=heading_log, def_col_width=12,
                               # background_color='black',
                               auto_size_columns=False,
                               display_row_numbers=True,
                               justification='center',
                               num_rows=20,
                               # alternating_row_color='black',
                               key='Table_Logs',
                               row_height=25, visible=True),
                       col_canvas,
                       sg.Col(columns_logs,key = 'Col_log',visible = True),
                       ]]
#+++++++++++++++++++++++++++++++++++++Thermal Reservoir tab ++++++++++++++++++++
TR_ROWS, TR_COLS, TR_HEADINGS = 9, 4, ('Min','Mid','Max','Method')
Demo_input_data = [[1000,1000,1000],[22000,22000,22000],[15,15,15],[2574,2623,2673],[778,824,876],
                   [4180,4180,4180],[31,36,43],[30,31,33],[180,520,620]]
Demo_method = ['Constant','Constant','Constant','Uniform','Normal','Laplace','Triangular','Triangular','Triangular']
# header_list = COL_HEADINGS
headings = list(TR_HEADINGS)
#name of rows, also used as keys for Combo
row_name = ['Water Density kg/m³','Area km²','Surface Temperature °C','Rock Density kg/m³','Specific Heat of Rock J/(kg·K)',
              'Specific Heat of Water J/(kg·K)','Reservoir Temperature °C','Porosity %','Thickness m']
TR_input_col = sg.Col(
    [[sg.Text(row_name[r],size = (12,2))]+[sg.Input(Demo_input_data[r][c], justification='c',size= (10,2), key=('TR'+str(r), 'TR'+str(c))) for c in range(TR_COLS-1)] +
    [sg.Combo(values=('Normal','Laplace','Triangular','Uniform','Constant'), default_value=Demo_method[r], key=row_name[r])] for r in range(TR_ROWS)]) #row names are  also used as keys

TR_output = sg.Frame('Geothermal Energy Result (J)',
                    [[sg.Text('                  Qtotal'),sg.Text('            Qrock'),sg.Text('            Qwater')],
                     [sg.Text('Median: \n Mean: \n std: \n Min: \n Max: \n 2.5%： \n 97.5%'),sg.Multiline(key='Q_total', size=(10,7)),sg.Multiline(key='Q_rock', size=(10,7)),sg.Multiline(key='Q_water', size=(10,7))]]
                    )
TR_left_col = sg.Col([[sg.Text('Geothermal Energy Estimating', font='Default 16')],
                    [sg.Text(' ' * 26, visible=True)] + [sg.Text(s, key=s, visible=True, enable_events=True, font='Default 12', size=(8, 1)) for s in TR_HEADINGS],
                    [TR_input_col],
                    [sg.Text('Iterations')]+[sg.Input(100000, key= 'TR_iter',justification='r',size= (10,1))]+[sg.Button('Run',key = 'Run_TR')],
                      [TR_output]
                      ])
figure_w, figure_h = 300, 300
TR_canvas = sg.Col([[sg.Canvas(size=(figure_w, figure_h), key='canvas_TR')],[sg.Canvas(key='controls_TR')]]) # cavas for TR
tab_TR_layout = [[TR_left_col,TR_canvas]]
#+++++++++++++++++++++++++++++++++++++Temperature Distribution 1-dimision tab ++++++++++++++++++++
columns_temp1D = sg.Col([[sg.Button('Import Layers 1D')],
                [sg.Text("Upper Boundary Temperature °C")],
                [sg.Input(0, key= 'upper_bd_val_1D',justification='c',size= (10,1))],
                [sg.Text("Bottom Boundary Condition")],
                [sg.Radio('Heat flow', "Botton_bound1D", default=True, size=(10,1), k='Heat_flow_bd_1D'), sg.Radio('Temperature', "Botton_bound1D", default=False, size=(10,1), k='temp_bd_1D')],
                [sg.Text("Boundary Condition Values in °C or W/m² ")],
                [sg.Input(0.05, key= 'bottom_bd_val_1D',justification='c',size= (10,1))],
                [sg.Button('Calc Temp 1D',disabled= True)],
                [sg.FileSaveAs('Export Temp 1D',key='Export_Temp_1D', initial_folder=initial_folder,enable_events = True,disabled=True,file_types = (('MS EXCEL', '.xlsx'),))],  # TODO: better names
                [sg.Button('Read and Run Multiple Profiles',disabled= False)],
            ]
                )

figure_w, figure_h = 400, 400
Temp1D_canvas = sg.Col([[sg.Canvas(size=(figure_w, figure_h), key='canvas_Temp1D')],[sg.Canvas(key='controls_Temp1D')]]) # cavas for temp1D
initial_layer = [[0,0,0]] #default values
heading_layer = ['Layer','Depth (m)','k (W/m·K)','A (W/m³)']
tab_temp1D_layout = [[sg.T('Temperature Distribution 1D',font='Default 18')],[sg.Table(values=initial_layer, headings=heading_layer, def_col_width=8,
                               # background_color='black',
                               auto_size_columns=False,
                               display_row_numbers=True,
                               justification='center',
                               num_rows=10,
                               # alternating_row_color='black',
                               key='Table_Temp1D',
                               row_height=25, visible=True),Temp1D_canvas,columns_temp1D]]

#+++++++++++++++++++++++++++++++++++++Temperature Distribution 2-dimision tab ++++++++++++++++++++
columns_temp2D = sg.Col([[sg.Button('Import Layers Data 2D')],
                [sg.Text("Upper Boundary Temperature °C")],
                [sg.Input(0, key= 'upper_bd_val_2D',justification='c',size= (10,1))],
                [sg.Text("Bottom Boundary Condition")],
                [sg.Radio('Heat flow', "Botton_bound2D", default=True, size=(10,1), k='Heat_flow_bd_2D'), sg.Radio('Temperature', "Botton_bound2D", default=False, size=(10,1), k='temp_bd_2D')],
                [sg.Button('Calc Temp 2D',disabled= True)],
                [sg.FileSaveAs('Export Temp 2D',key='Export_Temp_2D', initial_folder=initial_folder,enable_events = True,disabled=True,file_types = (('Comma-Separated Values', '.csv'),))]  # TODO: better names
            ])

figure_w, figure_h = 400, 400
Temp2D_canvas = sg.Col([[sg.Canvas(size=(figure_w, figure_h), key='canvas_Temp2D')],[sg.Canvas(key='controls_Temp2D')]]) # cavas for temp1D
initial_layer_2D = [[0,0]] #default values
heading_layer_2D = ['k (W/m·K)','A (W/m³)']
tab_temp2D_layout = [[sg.T('Temperature Distribution 2D',font='Default 18')],[sg.Table(values=initial_layer_2D, headings=heading_layer_2D, def_col_width=8,
                               # background_color='black',
                               auto_size_columns=False,
                               display_row_numbers=True,
                               justification='center',
                               num_rows=10,
                               # alternating_row_color='black',
                               key='Table_Temp2D',
                               row_height=25, visible=True),Temp2D_canvas,columns_temp2D]]

#+++++++++++++++++++++++++++++++++++++Tectono-Thermo Evolution single ++++++++++++++++++++
columns_subsidence_thermo_console = sg.Col([[sg.Button('Import Strata Data')],
                [sg.Button('Calc Subsidence',disabled= True)],
                [sg.Button('Calc Thermo-Evolution',disabled= True)],
                [sg.FileSaveAs('Export Subsidence',key='Export_Subsidence', initial_folder=initial_folder,enable_events = True,disabled=True,file_types = (('Comma-Separated Values', '.csv'),))], # TODO: better names
                [sg.FileSaveAs('Export Thermo-Evolution', key='Export_Themo_Evolution', initial_folder=initial_folder,
                                        enable_events=True, disabled=True,
                                            file_types = (('MS EXCEL', '.xlsx'),)
                )]])

figure_w, figure_h = 400, 400
Subsidence_canvas = sg.Col([[sg.Canvas(size=(figure_w, figure_h), key='canvas_subsidence')],[sg.Canvas(key='controls_subsidence')]]) # cavas for temp1D
#subsidence tables
Subsidence_init = [[0,0,0,0,0,0,0,0,0]] #default values
Subsidence_init_heading = ['Formation','Age','Depth','Paleobathymetry','Paleosealevel','Sand','Mud','Carbon','Coal']

#subsidences parameters
Sub_para_row, Sub_para_col, Sub_COL_HEADINGS, Sub_Col_ROWNAMES = 3, 4, ('Sand','Mud','Carbon', 'Coal'),('Phi0 (0-1)','c m⁻¹','Density kg/m³')
subsidence_keys = [['phi0_sand','phi0_mud','phi0_carbon','phi0_coal'],['c_sand','c_mud','c_carbon','c_coal'],['rho_sand','rho_mud','rho_carbon','rho_coal']]
Default_phi_values = [[0.6,0.5,0.6,0.9],[0.217e-3,0.515e-3,0.22e-3,0.7e-3],[2800,2400,2720,1800]]
col_sub_para =  [[sg.T(Sub_Col_ROWNAMES[r], size=(6, 2))] + [sg.Input(Default_phi_values[r][c], justification='c',size= (10,2),
                key= subsidence_keys[r][c]) for c in range(Sub_para_col)] for r in range(Sub_para_row)]

#thermo-evo parameters
thermo_para_row, thermo_para_col, thermo_COL_HEADINGS, thermo_Col_ROWNAMES = 4, 6, ('Depth m','Density kg/m³','TC W/(m·K)', 'A W/m³','Cp J/(kg·K)','Alpha °C^-1'),\
                                                                             ('Upper crust','Lower crust','Lithosphere\nmantle','Asthenosphere\nmantle')
default_thermo_evo_para = [[22000,2750,2.6,1.2e-6,1000,3.28e-5],[33000,2950,3.0,0.3e-6,1000,3.28e-5],[125000,3300,3.4,0.0003e-6,1000,3.28e-5],[6600000,3340,10.0,0.0003e-6,1000,3.28e-5]]

# header_list = COL_HEADINGS
thermo_col_headings = list(thermo_COL_HEADINGS)
col_thermo_evo_para =  [[sg.T(thermo_Col_ROWNAMES[r], size=(6, 2))] +[sg.Input(default_thermo_evo_para[r][c], justification='c',size= (10,2),
                                          key=('thermo'+str(r),'thermo'+ str(c))) for c in range(thermo_para_col)] for r in range(thermo_para_row)]
#input simulated table
sub_input_frame = sg.Frame('Input',
                           [[sg.Table(values=Subsidence_init, headings=Subsidence_init_heading, def_col_width=6,
                               # background_color='black',
                               auto_size_columns=False,
                               display_row_numbers=True,
                               justification='center',
                               num_rows=10,
                               # alternating_row_color='black',
                               key='subsidence',
                               row_height=25, visible=True)],
                        [sg.T('Tectonic subsidence parameters',font='Default 14')],
                        [sg.Text(' ' * 15, key='subpara_', visible=True)] + [
                                sg.Text(s, key='col_sub_hd' + str(s), visible=True, enable_events=True,
                                        font='Courier 11', size=(8, 1)) for
                                     i, s in enumerate(Sub_COL_HEADINGS)],
                        [sg.Col(col_sub_para,key = 'sub_paras_col',visible = True)],
                        [sg.Text('Water Density kg/m³'),sg.Input(1000,key = 'water_density',justification='r',size= (8,1)),
                          sg.Text('  Mantle Density kg/m³'),sg.Input(3300,key = 'mantle_density',justification='r',size= (8,1))
                          ,sg.Text('  PHI'),sg.Input(1,key = 'PHI',justification='r',size= (8,1))],
                        [sg.T('\nThermal evolution parameters',font='Default 14')],
                         [sg.Text(' ' * 15, key='thermopara_', visible=True)]+ [sg.Text(s, key='col_thermo_hd'+str(s), visible=True, enable_events=True, font='Courier 11', size=(8, 2)) for
                             i, s in enumerate(thermo_col_headings)],
                        [sg.Col(col_thermo_evo_para,key = 'sub_paras_col',visible = True)],
                        [sg.Text('Surface temperature °C'),sg.Input(0,key = 'T0_thermo',justification='r',size= (8,1))
                          ,sg.Text('Mantle Temperature °C'),sg.Input(1333,key = 'Tm_thermo',justification='r',size= (8,1))]
                        ])

tab_subsidence_layout = [[sg.T('Tectono-Thermo Evolution',font='Default 18')],
                         [sub_input_frame,Subsidence_canvas,columns_subsidence_thermo_console]]

#+++++++++++++++++++++++++++++++++++++Tectono-Thermo Evolution Batch ++++++++++++++++++++
Thermo_evo_batch_conso =sg.Col([
                                [sg.Button('Thermo Evo Folder',disabled = False)],
                                [sg.Button('Subsidence_Batch',disabled = True)],
                                [sg.Button('Thermo_Batch',disabled = True)]]) # functions
#simulated tables
col_sub_para_batch =  [[sg.T(Sub_Col_ROWNAMES[r], size=(6, 2))] + [sg.Input(Default_phi_values[r][c], justification='c',size= (10,2),
                key='batch_'+ subsidence_keys[r][c]) for c in range(Sub_para_col)] for r in range(Sub_para_row)]
col_thermo_evo_para_batch =  [[sg.T(thermo_Col_ROWNAMES[r], size=(6, 2))] +[sg.Input(default_thermo_evo_para[r][c], justification='c',size= (10,2),
                                          key=('batch_thermo'+str(r),'batch_thermo'+ str(c))) for c in range(thermo_para_col)] for r in range(thermo_para_row)]
sub_input_frame = sg.Frame('Input Batch',
                        [[sg.T('Tectonic subsidence parameters',font='Default 14')],
                        [sg.Text(' ' * 15, key='subpara_batch', visible=True)] + [
                                sg.Text(s, key='col_sub_batch' + str(s), visible=True, enable_events=True,
                                        font='Courier 11', size=(8, 1)) for
                                     i, s in enumerate(Sub_COL_HEADINGS)],
                        [sg.Col(col_sub_para_batch,key = 'sub_paras_batch',visible = True)],
                        [sg.Text('Water Density kg/m³'),sg.Input(1000,key = 'batch_water_density',justification='r',size= (8,1)),
                          sg.Text('  Mantle Density kg/m³'),sg.Input(3300,key = 'batch_mantle_density',justification='r',size= (8,1))
                          ,sg.Text('  PHI'),sg.Input(1,key = 'batch_PHI',justification='r',size= (8,1))],
                        [sg.T('\nThermal evolution parameters',font='Default 14')],
                         [sg.Text(' ' * 15, key='thermopara_batch', visible=True)]+ [sg.Text(s, key='col_thermo_batch'+str(s), visible=True, enable_events=True, font='Courier 11', size=(8, 2)) for
                             i, s in enumerate(thermo_col_headings)],
                        [sg.Col(col_thermo_evo_para_batch,key = 'sub_paras_col_batch',visible = True)],
                        [sg.Text('Surface temperature °C'),sg.Input(0,key = 'T0_thermo_batch',justification='r',size= (8,1))
                          ,sg.Text('Mantle Temperature °C'),sg.Input(1333,key = 'Tm_thermo_batch',justification='r',size= (8,1))]
                        ])
input_output_info_frame =sg.Frame('Input and output files',[[sg.Multiline(size=(50,20), font='Courier 8',key='input_dir',disabled = False)]+
                                                            [sg.Multiline(size=(50,20), font='Courier 8',key='output_dir',disabled = False)]])
tab_tectono_thermo_batch_layout = [[sg.T('Tectono-Thermo Evolution Batch',font='Default 18')],
                                   [input_output_info_frame],
                                   [sub_input_frame]+[Thermo_evo_batch_conso]]
#++++++++++++++++++++++++++++++++++++++Welcome tab+++++++++++++++++++++++++++++++

tab_welcome_layout = [[sg.T('Many thanks for using my code :)',font='Default 30')],
                      [sg.T('Let me know, when you find bugs.',font='Default 25')],
                      [sg.T('Written by Dr.Jie HU, Email: hujie@cdut.edu.cn',font='Default 25')],
                        [sg.T('Version 1.0',font='Default 25')]
                      ]

#++++++++++++++++++++++++++++++++++++++All layout+++++++++++++++++++++++++++++++
layout = [[sg.Text('Geothermal Tools', size=(50, 1), justification='center', font=("Helvetica", 20),
                   relief=sg.RELIEF_RIDGE, k='-TEXT HEADING-', enable_events=True)]]
layout += [[sg.TabGroup([[sg.Tab('Welcome',tab_welcome_layout),sg.Tab('TC', tab_TC_layout),sg.Tab('Well Log',tab_Well_log_layout),sg.Tab('Thermal Reservoir',tab_TR_layout),sg.Tab('Temperature 1D',tab_temp1D_layout),sg.Tab('Temperature 2D',tab_temp2D_layout)
                          ,sg.Tab('Tectono-Thermo Evo',tab_subsidence_layout),sg.Tab('Tectono-Thermo Evo Batch',tab_tectono_thermo_batch_layout)]], key='-TAB GROUP-')]]
# Create the window
window = sg.Window('Geothermal Tools', layout, default_element_size=(20, 1), element_padding=(1, 1), return_keyboard_events= True, resizable=False, finalize=True)
current_cell = (0, 0)
#===================================== Gui ends================================
#=====================================Main functions===========================
import_TC = False # the import_TC default to False
#
while True:
    event, values = window.read()
    print(event, values)
    # print("event:", event, "values: ", values)

#+++++++++++++++++++++++++++TC start+++++++++++++++++++++++++++++
    if event =='Calc TC':
        if import_TC== False:
            table = [[values[(row, col)] for col in range(MAX_COLS)] for row in range(MAX_ROWS)]
            tablepd = pd.DataFrame(table, columns=COL_HEADINGS, dtype=float)
        else:
            tablepd = impt_TC
            tablepd.columns = COL_HEADINGS
            #tablepd = pd.DataFrame(table, columns=COL_HEADINGS, dtype=np.float)

        # sg.popup_scrolled('your_table = [ ', ',\n'.join([str(table[i]) for i in range(MAX_ROWS)]) + '  ]', title='Copy your data from here', font='fixedsys', keep_on_top=True)
        #tablepd = pd.DataFrame(table,columns = COL_HEADINGS,dtype=np.float)
        T0 = values['T0']
        Unmixing_law = values['Unmixing_Law_Combo']
        Mixing_law = values['Mixing_Law_Combo']
        Fluid = values['Fluid_Combo']
        TC_s,TC_c,HF = corr_TC_HF(tablepd,T0,Unmixing_law,Mixing_law,Fluid)
        # tablepd['TCs\n(W/mK)'] = np.round(TC_s,2)
        # tablepd['TCc\n(W/mK)'] = np.round(TC_c,2)
        # tablepd['HF\n(mW/m2)'] = np.round(HF,2)

        tablepd.iloc[:,5] = np.round(TC_s,2)
        tablepd.iloc[:,6]= np.round(TC_c,2)
        tablepd.iloc[:,7]= np.round(HF, 2)
        # tablepd['HF\n(mW/m2)'] = np.round(HF,2)
        #if we havenot  import the TC data, update the table
        if import_TC== False:
            for i in range(MAX_ROWS):
                for j in range(MAX_COLS):
                    window[(i, j)].update(tablepd.iloc[i][j])
        else:
            window['Table_TC'].update(values=impt_TC[:].values.tolist()) # update the data in DataFrame
        window['Export_TC'].Update(disabled=False)  # enable the Export TC
    elif event == 'Import TC':
        try:  # avoid fail
            import_TC = True
            filename_TC = sg.popup_get_file(
                'filename to open', no_window=True, file_types=(("MS Excel", "*.xlsx"),))  # open file
            # print(filename_TC)
            impt_TC = pd.read_excel(filename_TC)     #read file
            # print(impt_TC)
            window['header_TC'].Update(visible=False)  # hide the simulated table
            window['Col0'].Update(visible=False)
            window['Table_TC'].Update(visible=True)
            for k in COL_HEADINGS:   # hide the simulated table header
                window[k].Update(visible=False)
        except:
            print('Please close the excel file')
            pass

    elif event =='Export_TC':
        print('Save')
        file_name_excel = values['Export_TC']
        try :
            tablepd.to_excel(file_name_excel,index = False)
        except:
            pass
#++++++++++++++++++++++++++++++++++++TC ends+++++++++++++++++++++++

    elif event == 'Import Logs':
        filename_Logs = sg.popup_get_file(
            'filename to open', no_window=True, file_types=(("MS Excel", "*.xlsx"),))  # open file
        # impt_Logs = pd.read_excel(filename_Logs)
        try:                # avoid fail
            print(filename_Logs)
            impt_Logs = pd.read_excel(filename_Logs)     #read file logs
            print(impt_Logs)
            # print('aa')
            window['Table_Logs'].update(values=impt_Logs[:].values.tolist())
            # window['Export Logs'].Update(disabled=False)  #
              #
            window['Plot Logs'].Update(disabled=False)  #
        except:
            pass

    elif event == 'Plot Logs':
        #-------plot code
        x = impt_Logs['Temperature']
        y = impt_Logs['Depth']
        fig0, axes = plt.subplots(1, 2,sharey=True,figsize=(5,5))
        axes[0].clear()
        axes[0].plot(x, y)
        axes[0].invert_yaxis()
        axes[0].set_title('Temperature')
        axes[0].set_xlabel('Temperature °C')
        axes[0].set_ylabel('Depth (m)')
        #---calc Gradient
        n_of_stp_logs = int(values['number_of_stp_logs'])
        x_for_grad = np.linspace(x[0],x[len(x)-1],n_of_stp_logs) #50 node
        from scipy import interpolate
        temp_func = interpolate.interp1d(x, y)
        y_for_grad = temp_func(x_for_grad)
        Grad = (x_for_grad[1:]-x_for_grad[:-1])/(y_for_grad[1:]-y_for_grad[:-1])*1000
        axes[1].step(np.array(Grad),np.array(y_for_grad[1:]),where = 'pre')
        axes[1].set_title('Gradient')
        axes[1].set_xlabel('Gradient °C/km')
        plt.subplots_adjust(wspace=0)
        plt.tight_layout()
        #---------put the figure into the canvas
        draw_figure_w_toolbar(window['canvas_logs'].TKCanvas, fig0, window['controls_cv'].TKCanvas)
        window['Fit Logs'].Update(disabled=False)
        window['Fit all'].Update(disabled=False)
    elif event == 'Fit Logs':
        number_segments = values['Number_of_segments']
        if number_segments ==1:  # if only 1 segment, use traditional fit
            z = np.round(np.polyfit(y, x, 1), 5) # fit results
            p = np.poly1d(z)     #fit function
            yHat = y
            xHat = p(yHat)
            breaks = z[1]
            slopes = z[0]
            # from sklearn.metrics import r2_score
            # rsq = np.round(r2_score(x, xHat), 4)
            _, _, r_value, _, _ = stats.linregress(y, x)
            rsq = round(r_value**2,4)
        else:
            import pwlf
            my_pwlf = pwlf.PiecewiseLinFit(y, x)
            #

            breaks = np.around(my_pwlf.fit(number_segments),1)
            slopes = np.around(my_pwlf.calc_slopes(),5)
            # predict for the determined points
            yHat = np.linspace(min(y), max(y), num=1000)
            xHat = my_pwlf.predict(yHat)
            rsq = np.around(my_pwlf.r_squared(),4) # r square Calculate the prediction variance
        axes[0].clear()
        axes[0].plot(xHat, yHat,label = 'Fitted')
        axes[0].scatter(x, y,s=3,c='r',label = 'Observed')
        axes[0].invert_yaxis()
        axes[0].set_title('Temperature')
        axes[0].set_xlabel('Temperature °C')
        axes[0].set_ylabel('Depth (m)')
        axes[0].legend()# replot the gradient
        plt.subplots_adjust(wspace=0)
        plt.tight_layout()
        #---------put the figure into the canvas
        draw_figure_w_toolbar(window['canvas_logs'].TKCanvas, fig0, window['controls_cv'].TKCanvas)
        #add information to fit output box
        mline: sg.Multiline = window['Fit_information']
        mline.update('slope: \n', append=False) # False is used to clear the previous text
        mline.update(slopes, append=True)
        mline.update('\nBreaks: \n', append=True)
        mline.update(breaks, append=True)
        mline.update('\nR square: \n', append=True)
        mline.update(rsq, append=True)
    elif event == 'Fit all':
        # fitting from 1 to 4
        all_rsq = []
        axes[0].clear()
        axes[0].scatter(x, y, s=3, c='r', label='Observed')
        axes[0].invert_yaxis()
        axes[0].set_title('Temperature')
        axes[0].set_xlabel('Temperature °C')
        axes[0].set_ylabel('Depth (m)')
        mline: sg.Multiline = window['Fit_information']
        mline.update('R square: \n', append=False) # False is used to clear the previous text
        for ns in (np.linspace(1,4,4,dtype = int)):
            if ns == 1:  # if only 1 segment, use traditional fit
                z = np.round(np.polyfit(y, x, 1), 5)  # fit results
                p = np.poly1d(z)  # fit function
                yHat = y
                xHat = p(yHat)
                breaks = z[1]
                slopes = z[0]
                # from sklearn.metrics import r2_score
                # rsq = np.round(r2_score(x, xHat), 4)
                _, _, r_value, _, _ = stats.linregress(y, x)
                rsq = round(r_value ** 2, 4)
            else:
                import pwlf
                my_pwlf = pwlf.PiecewiseLinFit(y, x)
                #
                breaks = np.around(my_pwlf.fit(ns), 1)
                slopes = np.around(my_pwlf.calc_slopes(), 5)
                # predict for the determined points
                yHat = np.linspace(min(y), max(y), num=1000)
                xHat = my_pwlf.predict(yHat)
                rsq = np.around(my_pwlf.r_squared(), 4)  # r square Calculate the prediction variance
            axes[0].plot(xHat, yHat, label=str(ns)+' Segments')
            all_rsq.append(rsq)
            axes[0].legend()  # replot the gradient
            plt.subplots_adjust(wspace=0)
            plt.tight_layout()
            # ---------put the figure into the canvas
            draw_figure_w_toolbar(window['canvas_logs'].TKCanvas, fig0, window['controls_cv'].TKCanvas)
            # add information to fit output box
            mline.update('\n'+str(ns)+' Segments:', append=True)
            # print(str(ns))
            mline.update(rsq, append=True)
    #
    elif event == 'Run_TR':
        All_TR_enter = [[values[('TR'+str(row), 'TR'+str(col))] for col in range(TR_COLS-1)] for row in range(TR_ROWS)]
        All_TR_values = pd.DataFrame(All_TR_enter,dtype = float)
        ALL_TR_method = [values[med] for med in row_name]

        Min_arr = All_TR_values[0] * 0.999  # to avoid the error
        Mean_arr = All_TR_values[1]
        Max_arr = All_TR_values[2] * 1.001
        iters = int(values['TR_iter'])

        all_rand_values = []
        # generate randon values
        for i in range(9):
            # gen_rand_values(df.loc[i][0],df.loc[i][1],df.loc[i][2],1000)
            all_rand_values.append(gen_rand_values(Min_arr[i], Mean_arr[i], Max_arr[i], ALL_TR_method[i], iters))
        # calc thermal energe
        TE_rand, Tr_rand, Tw_rand = thermal_energy(all_rand_values[0], all_rand_values[1], all_rand_values[2], \
                                                       all_rand_values[3], all_rand_values[4], all_rand_values[5], \
                                                       all_rand_values[6], all_rand_values[7], all_rand_values[8])
        # plot
        from matplotlib.figure import Figure
        # fig_TR = Figure(figsize=(10, 13), dpi=65)
        fig_TR= plt.figure(figsize=(10, 13),dpi=50)
        axes = fig_TR.subplots(4, 3)
        # print(axes)
        for i in range(9):
            k = i // 3  # the vertial position of ax
            j = i % 3  # the laternal position of ax
            plotter(axes[k][j], row_name[i], all_rand_values[i], Min_arr[i], Mean_arr[i], Max_arr[i], ALL_TR_method[i])

        # remove the useless plot
        axes[3][0].remove()
        axes[3][1].remove()
        axes[3][2].remove()
        axe_hist = fig_TR.add_subplot(4, 2, 7)
        axe_cdf = fig_TR.add_subplot(4, 2, 8)
        #
        counts, bins = np.histogram(TE_rand, bins=50)
        # # print(counts,bins)
        axe_hist.hist(bins[:-1], bins, weights=counts / iters)
        axe_hist.set_title('Geothermal resources histogram')
        axe_hist.set_xlabel('Values (J)')
        axe_hist.set_ylabel('Probability')
        from matplotlib.ticker import PercentFormatter

        axe_hist.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        # fill 95% confidence intervel
        TE_rand_sort = np.sort(TE_rand)
        TE_rand_first_2_5 = TE_rand_sort[int(0.025 * iters)]
        TE_rand_last_2_5 = TE_rand_sort[int(0.975 * iters)]
        axe_hist.axvspan(TE_rand_first_2_5, TE_rand_last_2_5, alpha=0.2, color='red')
        counts_lim = counts.max() / iters  # the max y
        axe_hist.set_ylim([0, counts_lim * 1.15])
        axe_hist.axhline(counts_lim * 1.05, color='k')
        axe_hist.text(TE_rand_first_2_5 / 2 + TE_rand_last_2_5 / 2, counts_lim * 1.06, '95%', ha='center')
        axe_hist.text(TE_rand_first_2_5 / 2 + TE_rand_sort[0] / 2, counts_lim * 1.06, '2.5%', ha='center')
        axe_hist.text(TE_rand_last_2_5 / 2 + TE_rand_sort[-1] / 2, counts_lim * 1.06, '2.5%', ha='center')

        axe_cdf.hist(bins[:-1], bins, weights=counts / iters, cumulative=True)
        # axe_cdf.plot(x,y)
        axe_cdf.set_title('Geothermal resources CDF')
        axe_cdf.set_xlabel('Values (J)')
        # axe_cdf.set_ylabel('Probability')
        axe_cdf.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        axe_cdf.axvspan(TE_rand_first_2_5, TE_rand_last_2_5, alpha=0.2, color='red')
        axe_cdf.set_ylim([0, 1.15])
        axe_cdf.axhline(1.05, color='k')
        axe_cdf.text(TE_rand_first_2_5 / 2 + TE_rand_last_2_5 / 2, 1.06, '95%', ha='center')
        axe_cdf.text(TE_rand_first_2_5 / 2 + TE_rand_sort[0] / 2, 1.06, '2.5%', ha='center')
        axe_cdf.text(TE_rand_last_2_5 / 2 + TE_rand_sort[-1] / 2, 1.06, '2.5%', ha='center')
        plt.subplots_adjust(top=0.95, bottom=0.05, left=0.08, right=0.98 )

        draw_figure_w_toolbar(window['canvas_TR'].TKCanvas, fig_TR, window['controls_TR'].TKCanvas)

        # output text
        # two decimals
        Out_med = np.median(TE_rand)
        oor = int(np.log10(Out_med))
        Out_med = int(str(int(Out_med) / 10 ** (oor - 3))[:3]) * 10 ** (
                    oor - 2) * 1.0  # use slice of str to get first 3 numbers
        # Out_med = round(Out_med/10**oor,3)*10**oor
        Out_mean = TE_rand.mean()
        oor = int(np.log10(Out_mean))
        Out_mean = int(str(int(Out_mean) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_std = TE_rand.std()
        oor = int(np.log10(Out_std))
        Out_std = int(str(int(Out_std) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_max = TE_rand.max()
        oor = int(np.log10(Out_max))
        Out_max = int(str(int(Out_max) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_min = TE_rand.min()
        oor = int(np.log10(Out_min))
        Out_min = int(str(int(Out_min) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        oor = int(np.log10(TE_rand_first_2_5))
        TE_rand_first_2_5 = int(str(int(TE_rand_first_2_5) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        oor = int(np.log10(TE_rand_last_2_5))
        TE_rand_last_2_5 = int(str(int(TE_rand_last_2_5) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        # put the values into output box
        mline: sg.Multiline = window['Q_total']
        mline.update(Out_med, append=False) # False is used to clear the previous text
        mline.update('\n', append=True)
        mline.update(Out_mean, append=True)
        mline.update('\n', append=True)
        mline.update(Out_std, append=True)
        mline.update('\n', append=True)
        mline.update(Out_min, append=True)
        mline.update('\n', append=True)
        mline.update(Out_max, append=True)
        mline.update('\n', append=True)
        mline.update(TE_rand_first_2_5, append=True)
        mline.update('\n', append=True)
        mline.update(TE_rand_last_2_5, append=True)

        # two decimals for rock
        Out_rock_med = np.median(Tr_rand)
        oor = int(np.log10(Out_rock_med))
        Out_rock_med = int(str(int(Out_rock_med) / 10 ** (oor - 3))[:3]) * 10 ** (
                    oor - 2) * 1.0  # use slice of str to get first 3 numbers
        # Out_rock_med = round(Out_rock_med/10**oor,3)*10**oor
        Out_rock_mean = Tr_rand.mean()
        oor = int(np.log10(Out_rock_mean))
        Out_rock_mean = int(str(int(Out_rock_mean) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_rock_std = Tr_rand.std()
        oor = int(np.log10(Out_rock_std))
        Out_rock_std = int(str(int(Out_rock_std) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_rock_max = Tr_rand.max()
        oor = int(np.log10(Out_rock_max))
        Out_rock_max = int(str(int(Out_rock_max) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_rock_min = Tr_rand.min()
        oor = int(np.log10(Out_rock_min))
        Out_rock_min = int(str(int(Out_rock_min) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Tr_rand = np.sort(Tr_rand)
        Tr_rand_first_2_5 = Tr_rand[int(len(Tr_rand) * 0.025)]
        Tr_rand_last_2_5 = Tr_rand[int(len(Tr_rand) * 0.975)]
        oor = int(np.log10(Tr_rand_first_2_5))
        Tr_rand_first_2_5 = int(str(int(Tr_rand_first_2_5) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        oor = int(np.log10(TE_rand_last_2_5))
        Tr_rand_last_2_5 = int(str(int(Tr_rand_last_2_5) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0

        # put the values into output box
        mline: sg.Multiline = window['Q_rock']
        mline.update(Out_rock_med, append=False) # False is used to clear the previous text
        mline.update('\n', append=True)
        mline.update(Out_rock_mean, append=True)
        mline.update('\n', append=True)
        mline.update(Out_rock_std, append=True)
        mline.update('\n', append=True)
        mline.update(Out_rock_min, append=True)
        mline.update('\n', append=True)
        mline.update(Out_rock_max, append=True)
        mline.update('\n', append=True)
        mline.update(Tr_rand_first_2_5, append=True)
        mline.update('\n', append=True)
        mline.update(Tr_rand_last_2_5, append=True)

        Out_water_med = np.median(Tw_rand)
        oor = int(np.log10(Out_water_med))
        Out_water_med = int(str(int(Out_water_med) / 10 ** (oor - 3))[:3]) * 10 ** (
                    oor - 2) * 1.0  # use slice of str to get first 3 numbers
        # two decimals
        Out_water_mean = Tw_rand.mean()
        oor = int(np.log10(Out_water_mean))
        Out_water_mean = int(str(int(Out_water_mean) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_water_std = Tw_rand.std()
        oor = int(np.log10(Out_water_std))
        Out_water_std = int(str(int(Out_water_std) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_water_max = Tw_rand.max()
        oor = int(np.log10(Out_water_max))
        Out_water_max = int(str(int(Out_water_max) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Out_water_min = Tw_rand.min()
        oor = int(np.log10(Out_water_min))
        Out_water_min = int(str(int(Out_water_min) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        Tw_rand = np.sort(Tw_rand)
        Tw_rand_first_2_5 = Tw_rand[int(len(Tw_rand) * 0.025)]
        Tw_rand_last_2_5 = Tw_rand[int(len(Tw_rand) * 0.975)]
        oor = int(np.log10(Tw_rand_first_2_5))
        Tw_rand_first_2_5 = int(str(int(Tw_rand_first_2_5) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0
        oor = int(np.log10(TE_rand_last_2_5))
        Tw_rand_last_2_5 = int(str(int(Tw_rand_last_2_5) / 10 ** (oor - 3))[:3]) * 10 ** (oor - 2) * 1.0

        # put the values into output box
        mline: sg.Multiline = window['Q_water']
        mline.update(Out_water_med, append=False)  # False is used to clear the previous text
        mline.update('\n', append=True)
        mline.update(Out_water_mean, append=True)
        mline.update('\n', append=True)
        mline.update(Out_water_std, append=True)
        mline.update('\n', append=True)
        mline.update(Out_water_min, append=True)
        mline.update('\n', append=True)
        mline.update(Out_water_max, append=True)
        mline.update('\n', append=True)
        mline.update(Tw_rand_first_2_5, append=True)
        mline.update('\n', append=True)
        mline.update(Tw_rand_last_2_5, append=True)

    elif event == 'Import Layers 1D':  # import layers for 1D
        try:
            filename_Temp1D = sg.popup_get_file(
                'filename to open', no_window=True, file_types=(("MS Excel", "*.xlsx"),))  # open file
            impt_layers1D = pd.read_excel(filename_Temp1D)  # read file logs
            window['Table_Temp1D'].update(values=impt_layers1D[:].values.tolist())
            # window['Export Logs'].Update(disabled=False)  #
            window['Calc Temp 1D'].Update(disabled=False)  #
        except:
            pass

    elif event == 'Calc Temp 1D':
        from Temp_func_1D import func_1D
        Layer_layer1D = impt_layers1D['Layer']
        Depth_layer1D = impt_layers1D['Depth']
        K_layer1D = impt_layers1D['K']
        A_layer1D = impt_layers1D['A']
        if values['Heat_flow_bd_1D'] == True:
            boundary_cond1D = 'HF'
        else:
            boundary_cond1D = 'Temp'
        upper_bound_1D = float(values['upper_bd_val_1D'])
        bottom_bound_1D = float(values['bottom_bd_val_1D'])
        xn1D,Tn1D = func_1D(Depth_layer1D,K_layer1D,A_layer1D,upper_bound_1D,bottom_bound_1D,boundary_cond1D)
        Temp_1D_prof = pd.DataFrame([])  # put values to DataFrame
        Temp_1D_prof['Depth'] = xn1D
        Temp_1D_prof['Temp'] = Tn1D
        # plot
        (_, file_name_1d_T) = os.path.split(filename_Temp1D)  # remove path
        figure_name_1DT = file_name_1d_T[:-5] # remove the extensional name
        fig_temp1D, ax_temp1D = plt.subplots(1, figsize=(4, 5))
        ax_temp1D.clear()
        ax_temp1D.plot(Tn1D, xn1D)
        ax_temp1D.invert_yaxis()
        ax_temp1D.set_title(figure_name_1DT)
        ax_temp1D.set_xlabel('Temperature °C')
        ax_temp1D.set_ylabel('Depth (m)')
        plt.tight_layout()
        # plt.show(bbox_inches = 'tight')

        # put in canvas
        draw_figure_w_toolbar(window['canvas_Temp1D'].TKCanvas, fig_temp1D, window['controls_Temp1D'].TKCanvas)
        window['Export_Temp_1D'].Update(disabled=False)

    elif event == 'Read and Run Multiple Profiles':
        fig_temp1D, ax_temp1D = plt.subplots(1, figsize=(4, 5))
        ax_temp1D.clear()
        ax_temp1D.invert_yaxis()
        ax_temp1D.set_xlabel('Temperature °C')
        ax_temp1D.set_ylabel('Depth (km)')
        plt.tight_layout()
        from Temp_func_1D import func_1D
        try:
            filename_Temp1D_multip = sg.popup_get_file(
                'filename to open', no_window=True, file_types=(("MS Excel", "*.xlsx"),))  # open file
            df_all_1D = pd.read_excel(filename_Temp1D_multip, sheet_name=None)  # read all sheets
        except:
            pass
        profile_names = list(df_all_1D.keys())
        (path_1D_multip, name_1D_multip) = os.path.split(filename_Temp1D_multip)# splite the path and name
        writer_fdez = pd.ExcelWriter(path_1D_multip+'/'+name_1D_multip[:-5]+'_Output.xlsx') # to save the results
        # print(path_1D_multip)
        for tp_name in profile_names:
            #read boundary conduction
            temp_profile = df_all_1D[tp_name]
            temp_boundary = temp_profile.head(1)
            temp_q0 = temp_boundary['q0'][0]
            temp_T0 = temp_boundary['T0'][0]
            #read layer
            temp_layers = pd.read_excel(filename_Temp1D_multip, sheet_name=tp_name,header = 2)
            temp_Depth_layer1D = np.array(temp_layers['Depth'])
            temp_K_layer1D = np.array(temp_layers['K'])
            temp_A_layer1D = np.array(temp_layers['A'])
            #thickness of each layer
            temp_thickness = temp_Depth_layer1D[1:]-temp_Depth_layer1D[:-1]
            temp_thickness = np.insert(temp_thickness,0,temp_Depth_layer1D[0])# insert the thickness of first layer
            #get the boundary conduction
            temp_qbbc = temp_q0 - np.sum(temp_thickness*temp_A_layer1D) #bottom boundary conduction

            temp_xn1D, temp_Tn1D = func_1D(temp_Depth_layer1D, temp_K_layer1D, temp_A_layer1D, temp_T0, temp_qbbc, 'HF')
            temp_1D_prof = pd.DataFrame([])  # put values to DataFrame
            temp_1D_prof['Depth'] = temp_xn1D
            temp_1D_prof['Temp'] = temp_Tn1D
            temp_1D_prof.to_excel(writer_fdez,sheet_name=tp_name,index=False)
            ax_temp1D.plot(temp_Tn1D, temp_xn1D/1000,label = tp_name)
        # #mantle adiabat T = 1250+0.5z
        # ma_z = np.linspace(0,200*1000,100)
        # ma_T = 1250+0.5/1000*ma_z
        # ax_temp1D.plot(ma_T,ma_z/1000)
        writer_fdez.save()
        writer_fdez.close()
        ax_temp1D.legend()
        # plt.savefig(filename_Temp1D[:-4]+'.pdf')
        # put in canvas
        draw_figure_w_toolbar(window['canvas_Temp1D'].TKCanvas, fig_temp1D, window['controls_Temp1D'].TKCanvas)
        window['Export_Temp_1D'].Update(disabled=False)

    elif event == 'Export_Temp_1D':
        file_name_export_1D = values['Export_Temp_1D']
        print(file_name_export_1D)
        try :
            Temp_1D_prof.to_excel(file_name_export_1D,index = False)
        except:
            pass

    elif event == 'Import Layers Data 2D':
        try:
            filename_Temp2D = sg.popup_get_file(
                'filename to open', no_window=True, file_types=(("MS Excel", "*.xlsx"),))  # open file
            impt_layers2D_ly = pd.read_excel(filename_Temp2D,sheet_name=0)  # read file layers
            impt_layers2D_para = pd.read_excel(filename_Temp2D, sheet_name=1)  # read file layers
            impt_layers2D_bhf = pd.read_excel(filename_Temp2D, sheet_name=2)  # read file layers

            window['Table_Temp2D'].update(values=impt_layers2D_para[:].values.tolist())
            # window['Export Logs'].Update(disabled=False)  #
            window['Calc Temp 2D'].Update(disabled=False)  #
        except:
            pass

    elif event == 'Calc Temp 2D':
        number_of_layers = int(len(impt_layers2D_ly.columns) / 2)
        #  from top to bottom
        layer_funcs = []
        for i in range(number_of_layers):
            x = impt_layers2D_ly.iloc[:, 2 * i].dropna()
            y = impt_layers2D_ly.iloc[:, 2 * i + 1].dropna()
            fc_temp = interpolate.interp1d(x, y)
            layer_funcs.append(fc_temp)

        #thermal parameters
        K_layer = impt_layers2D_para['K']  # thermal conductivity of each layer in W/mK
        A_layer = impt_layers2D_para['A']  # heat production in W/m3

        #boundary conditions
        upper_bounds = float(values['upper_bd_val_2D'])
        if values['Heat_flow_bd_2D'] == True:
            boundary_cond2D = 'HF'
        else:
            boundary_cond2D = 'Temp'
        # bottom_bd = pd.read_excel('Bottom_boundary_2D_test.xlsx')
        bottom_bounds = interpolate.interp1d(impt_layers2D_bhf.iloc[:, 0], impt_layers2D_bhf.iloc[:, 1])

        xlim=list(impt_layers2D_ly.iloc[:,0].dropna())[-1]# the last number in first column
        ylim = impt_layers2D_ly.iloc[:, -1][0]  # the first number in last column
        print(xlim,ylim)
        from Temp_func_2D import func_2D
        x_arr, y_arr, T_2D, hf = func_2D(xlim, ylim, layer_funcs, K_layer, A_layer, upper_bounds, bottom_bounds,
                                      boundary_cond2D)
        #plot
        fig_2d_temp = plt.figure(figsize=(10, 10),dpi=60)

        rect1 = [0.14, 0.75, 0.77, 0.2]  # [upper, lower, width, heigth] figure area （from 0 to 1）
        rect2 = [0.14, 0.3, 0.77, 0.4]
        rect3 = [0.14, 0.05, 0.77, 0.2]

        ax_shf = plt.axes(rect1)
        ax_shf.plot(x_arr, hf * 1000)
        ax_shf.set_xlim(x_arr[0], x_arr[-1])
        ax_shf.set_title('Surface Heat Flow mW/m²')

        ax_temp_field = plt.axes(rect2)
        # plot the layers
        for fc in layer_funcs:
            ax_temp_field.plot(x_arr, fc(x_arr), c='k')
        X_pcolor,Y_pcolor = np.meshgrid(x_arr,y_arr)
        c = ax_temp_field.pcolor(X_pcolor, Y_pcolor, T_2D, cmap='jet',shading='auto')
        # ax_temp_field.contour(X_pcolor, Y_pcolor, T_2D)
        ax_temp_field.invert_yaxis()
        cbar_ax = fig_2d_temp.add_axes(
            [0.92, 0.3, 0.015, 0.4])  # x,y of bottom left point, width, heigth of the colorbar
        fig_2d_temp.colorbar(c, cax=cbar_ax)
        ax_temp_field.set_title('Temperature Distribution °C')
        ax_temp_field.set_ylabel('Depth m')
        ax_bhf = plt.axes(rect3)
        ax_bhf.plot(x_arr, bottom_bounds(x_arr) * 1000)
        ax_bhf.set_xlim(x_arr[0], x_arr[-1])
        ax_bhf.set_title('Bottom Heat Flow mW/m²')
        ax_bhf.set_xlabel('Distance m')

        # put in canvas
        draw_figure_w_toolbar(window['canvas_Temp2D'].TKCanvas, fig_2d_temp, window['controls_Temp2D'].TKCanvas)
        window['Export_Temp_2D'].Update(disabled=False)

    elif event =='Export_Temp_2D':
        file_name_export_2D = values['Export_Temp_2D']
        print(file_name_export_2D)
        # X_2d_save = X_pcolor.stack()
        X_2d_save = X_pcolor.flatten() # change to 1D
        Y_2d_save = Y_pcolor.flatten()  # change to 1D
        T_2D_save = T_2D.flatten()
        Temp_2D_save = pd.DataFrame([])
        Temp_2D_save['x'] = X_2d_save
        Temp_2D_save['y'] = Y_2d_save
        Temp_2D_save['T'] = T_2D_save
        try:
            Temp_2D_save.to_csv(file_name_export_2D, index=False)
        except:
            pass

    elif event == 'Import Strata Data':
        try:
            filename_subsidence = sg.popup_get_file(
                'filename to open', no_window=True, file_types=(("Comma-Separated Values", "*.csv"),))  # open file
            impt_strata = pd.read_csv(filename_subsidence)  # read file logs
            window['subsidence'].update(values=impt_strata[:].values.tolist()) # present the data
            # window['Export Logs'].Update(disabled=False)  #
            window['Calc Subsidence'].Update(disabled=False)  #
        except:
            pass

    elif event == 'Calc Subsidence':
        from subsidence_single import cal_subsidence
        subsidc = cal_subsidence(filename_subsidence,float(values['phi0_sand']),float(values['phi0_mud']), float(values['phi0_carbon']), float(values['phi0_coal']),
                                 float(values['c_sand']), float(values['c_mud']), float(values['c_carbon']), float(values['c_coal']), float(values['rho_sand']), float(values['rho_mud']),
                       float(values['rho_carbon']), float(values['rho_coal']), float(values['water_density']), float(values['mantle_density']),float(values['PHI']))
        (_,file_name_sub_no_path) = os.path.split(filename_subsidence) # remove path
        (file_no_extension,_) = os.path.splitext(file_name_sub_no_path) #remove Filename Extension
        fig_sub = plt.figure(figsize=(7,5),dpi=65)
        ax_sub = fig_sub.add_subplot(1, 1, 1)
        ax_sub.plot(subsidc['age'],subsidc['subsidence'] )
        ax_sub.invert_xaxis()
        ax_sub.invert_yaxis()
        #    ax.legend(('observed subsidence','Fitted subsidence'))
        ax_sub.set_title(file_no_extension)
        ax_sub.set_xlabel('Age (Ma)')
        ax_sub.set_ylabel('Tectonic Subsidence (m)')

        #put in canvas
        draw_figure_w_toolbar(window['canvas_subsidence'].TKCanvas, fig_sub, window['controls_subsidence'].TKCanvas)

        window['Calc Thermo-Evolution'].Update(disabled=False)  #enable the Thermo-Evolution botton
        window['Export_Subsidence'].Update(disabled=False)  # enable the Thermo-Evolution botton

    elif event == 'Export_Subsidence':
        file_name_export_subsidence = values['Export_Subsidence'] # get file name
        try:
            subsidc.to_csv(file_name_export_subsidence, index=False)
        except:
            pass

    elif event == 'Calc Thermo-Evolution':
        from Thermo_evo import thermo_ev
        table_thermo = [[values[('thermo' + str(row), 'thermo' + str(col))] for col in range(thermo_para_col)] for row in range(thermo_para_row)]
        table_thermo_pd = pd.DataFrame(table_thermo, columns=thermo_col_headings, dtype=float) #read default values
        # read values
        T0_thermo = float(values['T0_thermo'])
        Tm_thermo = float(values['Tm_thermo'])
        depth_thermo = table_thermo_pd['Depth m']
        density_thermo = table_thermo_pd['Density kg/m³']
        TC_thermo = table_thermo_pd['TC W/(m·K)']
        A_thermo = table_thermo_pd['A W/m³']
        Cp_thermo = table_thermo_pd['Cp J/(kg·K)']
        Alpha_thermo = table_thermo_pd['Alpha °C^-1']
        # ('Depth m', 'Density kg/m³', 'TC W/(m·K)', 'A W/m³', 'Cp J/(kg·K)', 'Alpha °C^-1'), \
        # ('Upper crust', 'Lower crust', 'Lithosphere\nmantle', 'Asthenosphere\nmantle')
        # sg.popup('Wait a moment, the program is running :)')  # Remind the user
        #calc strain rate
        # Y,T,age_new,beta,HF_for_use,strain_rate= thermo_ev(subsidc,T0_thermo, Tm_thermo, depth_thermo, density_thermo, TC_thermo, A_thermo, Cp_thermo,Alpha_thermo,
        #           float(values['water_density']))
        Y, T, age_new, beta, HF_for_use, strain_rate, op_age = thermo_ev(subsidc,T0_thermo, Tm_thermo, depth_thermo, density_thermo, TC_thermo, A_thermo, Cp_thermo,Alpha_thermo,
                  float(values['water_density']))
        #subsidence
        fig_thermo, axes_thermo = plt.subplots(2, 2, sharex=True, figsize=(10, 8),dpi=65)
        axes_thermo[0][0].plot(subsidc['age'],subsidc['subsidence'])
        axes_thermo[0][0].plot(age_new, Y)
        axes_thermo[0][0].legend(('Observed subsidence', 'Fitted subsidence'))
        # axes_thermo[0][0].set_title(file_no_extension)
        axes_thermo[0][0].set_title('Subsidence')
        axes_thermo[0][0].invert_xaxis()
        axes_thermo[0][0].invert_yaxis()
        axes_thermo[0][0].set_ylabel('Depth m')
        #Heat flow
        HF_pd = pd.DataFrame([])
        HF_pd['age'] = age_new
        HF_pd['heat flow'] = HF_for_use
        axes_thermo[0][1].plot(age_new, HF_for_use)
        axes_thermo[0][1].set_ylabel('Heat flow mW/m²')
        axes_thermo[0][1].set_title('Heat flow')
        # axes_thermo[0][1].invert_xaxis()
        # Strain rate
        strain_ratepd = pd.DataFrame([])
        strain_ratepd['age'] = op_age
        strain_ratepd['strain rate'] = strain_rate
        axes_thermo[1][0].plot(op_age, strain_rate)
        axes_thermo[1][0].set_ylabel('Strain rate lg(s^-1)')
        axes_thermo[1][0].set_title('Strain rate')
        # axes_thermo[1][0].invert_xaxis()

        # beta
        beta_pd = pd.DataFrame([])
        beta_pd['age'] = age_new
        beta_pd['beta'] = beta
        axes_thermo[1][1].plot(age_new, beta)
        axes_thermo[1][1].set_ylabel('Beta')
        axes_thermo[1][1].set_title('Stretching factor')
        # axes_thermo[1][1].invert_xaxis()

        plt.tight_layout()
        draw_figure_w_toolbar(window['canvas_subsidence'].TKCanvas, fig_thermo, window['controls_subsidence'].TKCanvas)
        window['Export_Themo_Evolution'].Update(disabled=False)  # enable the Thermo-Evolution botton

    elif event == 'Export_Themo_Evolution':
        file_name_export_thermo_evo = values['Export_Themo_Evolution']  # get file name
        Thermo_evo_writer = pd.ExcelWriter(file_name_export_thermo_evo)
        try:
            subsidc.to_excel(Thermo_evo_writer, sheet_name = 'Subsidence' , index=False)
            HF_pd.to_excel(Thermo_evo_writer, sheet_name = 'Heat flow' , index=False)
            strain_ratepd.to_excel(Thermo_evo_writer, sheet_name = 'Strain rate' , index=False)
            Thermo_evo_writer.close()
        except:
            pass

    # elif event == 'Thermo Evo Folder':
    #     print('hj')
    elif event == 'Thermo Evo Folder':
        try:
            # thermo_evo_folder = values['thermo_ev_folder_batch']  # get folder name, change to inputdir
            thermo_evo_folder = sg.popup_get_folder('filename to open', no_window=True)
            names_all = os.listdir(thermo_evo_folder)
            mline: sg.Multiline = window['input_dir']
            mline.update('', append=False)  # False is used to clear the previous text
            for pf in names_all:
                mline.update(pf, append=True)
                mline.update('\n', append=True)
            window['Subsidence_Batch'].Update(disabled=False)
        except:
            pass

    elif event == 'Subsidence_Batch':
        from subsidence_single import cal_subsidence
        #make output dir
        input_name_base_name = os.path.basename(thermo_evo_folder) #base name
        fzd_name = os.path.dirname(thermo_evo_folder) # path
        out_dir = fzd_name+'/'+ input_name_base_name+'_output'
        try:
            os.mkdir(out_dir)            # make output dir
        except:
            print('az')
            pass
        mline: sg.Multiline = window['output_dir']
        mline.update('Subsidence Finished Well \n', append=False)  # False is used to clear the previous text
        subsidence_all = []
        for iz in names_all:
            try:
                filename_subsidence_batch_temp = iz
                os.chdir(thermo_evo_folder)  # change to the input dir
                subsidc_temp = cal_subsidence(filename_subsidence_batch_temp, float(values['batch_phi0_sand']), float(values['batch_phi0_mud']),
                                         float(values['batch_phi0_carbon']), float(values['batch_phi0_coal']),
                                         float(values['batch_c_sand']), float(values['batch_c_mud']), float(values['batch_c_carbon']),
                                         float(values['batch_c_coal']), float(values['batch_rho_sand']), float(values['batch_rho_mud']),
                                         float(values['batch_rho_carbon']), float(values['batch_rho_coal']), float(values['batch_water_density']), float(values['mantle_density']),
                                         float(values['batch_PHI']))
                subsidence_all.append(subsidc_temp)
                mline.update(iz, append= True)  # False is used to clear the previous text
                mline.update('\n', append= True)  # False is used to clear the previous text
                os.chdir(out_dir) # changed to output dir, using absolute path
                title = iz[:-4]
                subsidc_temp.to_csv(iz[:-4]+'.csv', index=False)
                #figs to save
                fig_sub_bat = plt.figure(figsize=(8, 5))
                ax_sub_bat = fig_sub_bat.add_subplot(1, 1, 1)
                ax_sub_bat.plot(subsidc_temp['age'],subsidc_temp['subsidence'] )
                ax_sub_bat.invert_xaxis()
                ax_sub_bat.invert_yaxis()
                ax_sub_bat.set_xlabel('Age (Ma)')
                ax_sub_bat.set_ylabel('Subsidence (m)')
                #    ax.legend(('observed subsidence','Fitted subsidence'))
                ax_sub_bat.set_title(title)
                fig_name_sub_bat = title + '.pdf'
                plt.savefig(fig_name_sub_bat)
                plt.close()
            except:
                print('No'+iz)

        window['Thermo_Batch'].Update(disabled=False)
    elif event == 'Thermo_Batch':
        from Thermo_evo import thermo_ev

        table_thermo_batch = [[values[('batch_thermo' + str(row), 'batch_thermo' + str(col))] for col in range(thermo_para_col)] for row
                        in range(thermo_para_row)]
        table_thermo_pd_batch = pd.DataFrame(table_thermo_batch, columns=thermo_col_headings, dtype=float)  # read default values
        # read values
        T0_thermo_batch = float(values['T0_thermo_batch'])
        Tm_thermo_batch = float(values['Tm_thermo_batch'])
        depth_thermo_batch = table_thermo_pd_batch['Depth m']
        density_thermo_batch = table_thermo_pd_batch['Density kg/m³']
        TC_thermo_batch = table_thermo_pd_batch['TC W/(m·K)']
        A_thermo_batch = table_thermo_pd_batch['A W/m³']
        Cp_thermo_batch = table_thermo_pd_batch['Cp J/(kg·K)']
        Alpha_thermo_batch = table_thermo_pd_batch['Alpha °C^-1']
        # ('Depth m', 'Density kg/m³', 'TC W/(m·K)', 'A W/m³', 'Cp J/(kg·K)', 'Alpha °C^-1'), \
        # ('Upper crust', 'Lower crust', 'Lithosphere\nmantle', 'Asthenosphere\nmantle')
        os.chdir(out_dir)
        mline.update('Thermal Evolution Finished Well \n', append=False)  # False is used to clear the previous text
        for zze in range(len(subsidence_all)):
            Y, T, age_new, beta, HF_for_use, strain_rate, op_age = thermo_ev(subsidence_all[zze], T0_thermo_batch, Tm_thermo_batch, depth_thermo_batch,
                                                                     density_thermo_batch, TC_thermo_batch, A_thermo_batch, Cp_thermo_batch,
                                                                     Alpha_thermo_batch,
                                                                     float(values['batch_water_density']))
            # subsidence
            fig_thermo, axes_thermo = plt.subplots(4, 1, sharex=True, figsize=(5, 10), dpi=65)
            thermo_ev_title_tp = names_all[zze][:-4]
            axes_thermo[0].plot(subsidence_all[zze]['age'], subsidence_all[zze]['subsidence'])
            axes_thermo[0].plot(age_new, Y)
            subsidence_fitted_df = pd.DataFrame([])
            subsidence_fitted_df['age'] = age_new
            subsidence_fitted_df['subsidence'] = Y
            axes_thermo[0].legend(('Observed subsidence', 'Fitted subsidence'))
            axes_thermo[0].set_title(thermo_ev_title_tp)
            axes_thermo[0].invert_xaxis()
            axes_thermo[0].invert_yaxis()
            axes_thermo[0].set_ylabel('Depth m')
            # Heat flow
            HF_pd = pd.DataFrame([])
            HF_pd['age'] = age_new
            HF_pd['heat flow'] = HF_for_use
            axes_thermo[1].plot(age_new, HF_for_use)
            axes_thermo[1].set_ylabel('Heat flow mW/m²')
            # axes_thermo[1].set_title('Heat Flow')
            # Strain rate
            strain_ratepd = pd.DataFrame([])
            strain_ratepd['age'] = op_age
            strain_ratepd['strain rate'] = strain_rate
            axes_thermo[2].plot(op_age, strain_rate)
            axes_thermo[2].set_ylabel('Strain rate s^-1')
            # axes_thermo[2].set_title('Strain rate')

            # Beta
            betapd = pd.DataFrame([])
            betapd['age'] = age_new
            betapd['beta'] = beta
            axes_thermo[3].plot(age_new, beta)
            axes_thermo[3].set_ylabel('Stretching factor')
            # axes_thermo[3].set_title('Stretching factor')
            axes_thermo[3].set_xlabel('Age (Ma)')

            mline.update(names_all[zze], append=True)  # False is used to clear the previous text
            mline.update('\n', append=True)  # False is used to clear the previous text

            plt.tight_layout()
            fig_name_sub_bat = thermo_ev_title_tp + '.pdf'
            plt.savefig(fig_name_sub_bat)
            Thermo_evo_writer_tp = pd.ExcelWriter(thermo_ev_title_tp+'.xlsx')
            try:
                subsidence_all[zze].to_excel(Thermo_evo_writer_tp, sheet_name='Subsidence observed', index=False)
                subsidence_fitted_df.to_excel(Thermo_evo_writer_tp,sheet_name='Subsidence Fitted', index=False)
                HF_pd.to_excel(Thermo_evo_writer_tp, sheet_name='Heat flow', index=False)
                strain_ratepd.to_excel(Thermo_evo_writer_tp, sheet_name='Strain rate', index=False)
                betapd.to_excel(Thermo_evo_writer_tp, sheet_name='Beta', index=False)
                Thermo_evo_writer_tp.close()
            except:
                pass
            sg.one_line_progress_meter('Batch',zze+1,len(subsidence_all))

    elif event in ('Exit', None):
        break

window.close()