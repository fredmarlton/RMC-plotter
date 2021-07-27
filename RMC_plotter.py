

# A simple and effective plotting code for assessing the fits from RMCProfile refinements


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.widgets import MultiCursor
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os
import fnmatch

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = '11'


# From local python modules
# Taken from https://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
import cm2inch as c2i

def TOF_2_Q(TOF_data, Lpath, two_theta):

    h_bar = 6.626176e-34
    m_n = 1.67495e-27
    sintheta = np.sin((np.pi/180)*two_theta/2)
    lambda_tof = (10**4)*h_bar*TOF_data/(m_n*Lpath) 

    Q = 4*np.pi*sintheta/(lambda_tof)

    return Q

# Useful function for seeing values on dual axis plots
def make_format(current, other):
    # current and other are axes
    def format_coord(x, y):
        # x, y are data coordinates
        # convert to display coords
        display_coord = current.transData.transform((x,y))
        inv = other.transData.inverted()
        # convert back to data coords with respect to ax
        ax_coord = inv.transform(display_coord)
        coords = [ax_coord, (x, y)]
        return ('Left: {:<40}    Right: {:<}'
                .format(*['({:0.3e}, {:0.3e})'.format(x, y) for x,y in coords]))
    return format_coord

def no_space(a):
    return a.replace(" ","")
remove_space = np.vectorize(no_space)

def middle_text(axis, the_text):
    # build a rectangle in axes coords
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    return axis.text(0.5*(left+right), 0.5*(bottom+top), the_text,
    horizontalalignment='center',
    verticalalignment='center',
    fontsize=20, color='red',
    transform=axis.transAxes,
    multialignment='left')

Bragg_text='No Bragg file'

def rmcplot(path, stem, pdf_type, recip_type, **kwargs):

    """
    ----
    Main Parameters:
    path = RMC refinement folder \n
    stem = RMC files stem name \n
    pdf_type = "ft_fq", "gr" \n
    recip_type =  "fq", "sq" \n
    \n
    ----
    Example: \n
    path = r"C:\\RMCProfile\\RMCProfile_package\\tutorial\\ex_6\\rmc_neutron\\run" \n
    stem_name = "gapo4_neutron" \n
    pdf_type = "gr" \n
    recip_type = "sq" \n
    ----
    kwargs: \n
    partials="Combine_partials" # When specified requires user input of partials, which is useful for mixed site materials \n
    partials_matrix=partials # Each column of this matrix corresponds to a partial pair, with the first column being the x-axis \n
    partials_labels=labels # Labels for each column of the partials_matrix \n
    recip_diff_fact=float # If specified will multiply the recip. difference curve by float for magnification \n
    pdf_diff_fact=float # If specified will multiply the pdf difference curve by float for magnification \n
    """

    path_stem = path + "\\" + stem

    # fig = plt.figure(figsize=c2i.cm2inch(60,30)) #width, height
    fig = plt.figure(figsize=c2i.cm2inch(40,20)) #width, height
    ax = plt.subplot(221) # PDF 
    ax1 = plt.subplot(223, sharex=ax) # PDF partials
    ax2 = plt.subplot(324) # Bragg
    ax3 = plt.subplot(322) # Chi2 data
    ax2_chi = ax3.twinx() # Chi2 data
    ax4 = plt.subplot(326) # Recip data



 ## ----- PDF ----- 
    if pdf_type=="ft_fq":
        PDF_data_file = path_stem + "_FT_XFQ1.csv"
    elif pdf_type=="gr":
        PDF_data_file = path_stem + "_PDF1.csv"
    elif pdf_type=="gr_x":
        PDF_data_file = path_stem + "_Xray_PDF1.csv"

    ax.set_xlabel(r'r($\AA$)')
    ax.set_ylabel('G(r)')
    ax.set_title('PDF fit')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    if os.path.exists(PDF_data_file):
        pdf_data = np.loadtxt(PDF_data_file, delimiter=",", skiprows=1)
        pdf_labels = np.genfromtxt(PDF_data_file, dtype='str', delimiter=",", max_rows=1)

        # pdf_rfd = np.loadtxt(path_stem+'.gr', skiprows=2) #original raw data file

        pdf_labels = remove_space(pdf_labels)

        diff_zoom = 1.0; diff_label='diff'
        if kwargs.get("pdf_diff_fact"):
            diff_zoom = kwargs.get("pdf_diff_fact")
            diff_label='diff' + 'x' + str(diff_zoom)

        pdf_r    = pdf_data[:,0]
        pdf_fit  = pdf_data[:,1]
        pdf_data = pdf_data[:,2]
        pdf_diff = pdf_data - pdf_fit
        pdf_diff = pdf_diff*diff_zoom
        diff_offset = abs(min(min(pdf_fit),min(pdf_data))) + abs(max(pdf_diff[100:])) + 0.05
        pdf_diff = pdf_diff - diff_offset

        # Black line through difference
        ax.axhline(-diff_offset, color='k')

        # ax.plot(pdf_rfd[:,0], pdf_rfd[:,1], 'bo-', markerfacecolor='white', label='.gr', markersize=4)
        ax.plot(pdf_r,pdf_data, 'bo-', markerfacecolor='white', label=pdf_labels[2], markersize=4)
        ax.plot(pdf_r,pdf_fit,  color='red', label=pdf_labels[1])
        ax.plot(pdf_r,pdf_diff,  color='green', label=diff_label)
        ax.set_xlim(left = min(pdf_r)-0.1, right= max(pdf_r)+0.1)
        ax.set_ylim(bottom=min(pdf_diff))
        ax.legend(loc=2, ncol=3)

        if kwargs.get('plot_baseline'):
            num_density = kwargs.get("num_density")
            bl_r = -pdf_r*4*np.pi*num_density
            # ax.plot(pdf_r, bl_r, 'k', zorder=-1000)
            ax.plot(pdf_r, bl_r, 'k', zorder=1000)


    else:
        middle_text(ax, 'No PDF fit')


    
 ## ----- Recip -----
    if recip_type=="fq":
        Recip_Data_File = path_stem + "_FQ1.csv"
    elif recip_type=="sq":
        Recip_Data_File = path_stem + "_SQ1.csv"
    
    ax4.set_xlabel(r'Q($\AA^{-1}$)')
    ax4.set_ylabel('F(Q),S(Q),etc.')
    ax4.set_title('Reciprocal space fit')
    ax4.xaxis.set_minor_locator(AutoMinorLocator())
    ax4.yaxis.set_minor_locator(AutoMinorLocator())

    
    if os.path.exists(Recip_Data_File):
        rec_data = np.loadtxt(Recip_Data_File, delimiter=",", skiprows=1)
        rec_labels = np.genfromtxt(Recip_Data_File, dtype='str', delimiter=",", max_rows=1)

        rec_labels = remove_space(rec_labels) 

        diff_zoom = 1.0; diff_label='diff'
        if kwargs.get("recip_diff_fact"):
            diff_zoom = kwargs.get("recip_diff_fact")
            diff_label='diff' + 'x' + str(diff_zoom)

        rec_q    = rec_data[:,0]
        rec_fit  = rec_data[:,1]
        rec_data = rec_data[:,2]
        rec_diff = rec_data - rec_fit
        rec_diff = rec_diff*diff_zoom
        diff_offset = abs(min(min(rec_fit),min(rec_data))) + abs(max(rec_diff)) + 0.05
        rec_diff = rec_diff - diff_offset

        # Black line through difference
        ax4.axhline(-diff_offset, color='k')

        # ax4.plot(rec_q,rec_data, 'bo-', markerfacecolor='white', label=rec_labels[2], markersize=4)
        ax4.plot(rec_q,rec_data, 'b-', lw=1.5, label=rec_labels[2], markersize=4)
        ax4.plot(rec_q,rec_fit,  color='red', lw=1, label=rec_labels[1])
        ax4.plot(rec_q,rec_diff,  color='green', label=diff_label)
        ax4.set_xlim(min(rec_q)-0.1, max(rec_q)+0.1)
        ax4.legend(loc=1)
    else:
        middle_text(ax4, 'No Reciprocal space fit')



 ## ----- Bragg -----
    Bragg_Data_File = path_stem + "_bragg.csv"

    ax2.set_xlabel('2Theta, TOF, etc.')
    ax2.set_ylabel('Intensity')
    ax2.set_title('Bragg fit')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())

    if os.path.exists(Bragg_Data_File):
        # bragg_data = np.loadtxt(Bragg_Data_File, delimiter=",", skiprows=1)
        bragg_data = np.loadtxt(Bragg_Data_File, delimiter=",", skiprows=0)
        # bragg_labels = np.genfromtxt(Bragg_Data_File, dtype='str', delimiter=",", max_rows=1)

        # bragg_labels = remove_space(bragg_labels) 

        diff_zoom = 1.0; diff_label='diff'
        if kwargs.get("recip_diff_fact"):
            diff_zoom = kwargs.get("recip_diff_fact")
            diff_label='diff' + 'x' + str(diff_zoom)

        bragg_x    = bragg_data[:,0]
        bragg_exp  = bragg_data[:,1]
        bragg_fit  = bragg_data[:,2]
        bragg_diff = bragg_exp - bragg_fit
        bragg_diff = bragg_diff*diff_zoom
        max_bragg = max(max(bragg_exp),max(bragg_fit))
        min_bragg = min(min(bragg_exp),min(bragg_fit))
        diff_offset= - (min_bragg - max(bragg_diff) - 0.03*(max_bragg-min_bragg))
        bragg_diff = bragg_diff - diff_offset

        if kwargs.get('plot_Bragg_in_Q'):
            bragg_q = TOF_2_Q(bragg_x, kwargs.get('L_path'), kwargs.get('two_theta'))
            bragg_x = bragg_q
            ax2.set_xlabel(r'Q($\AA^{-1}$)')

        if kwargs.get('plot_Bragg_in_d'):
            bragg_q = TOF_2_Q(bragg_x, kwargs.get('L_path'), kwargs.get('two_theta'))
            bragg_x = 2*np.pi/bragg_q
            ax2.set_xlabel(r'd($\AA$)')

        # Black line through difference
        ax2.axhline(-diff_offset, color='k')

        # ax2.plot(bragg_x,bragg_exp, 'b', marker='x', markersize=4, label=bragg_labels[1], linewidth=1.5)
        ax2.plot(bragg_x,bragg_exp, 'b', marker='x', markersize=3, label='data', linewidth=1.5)
        # ax2.plot(bragg_x,bragg_fit,  color='red', label=bragg_labels[2])
        ax2.plot(bragg_x,bragg_fit,  color='red', label='fit')
        ax2.plot(bragg_x,bragg_diff,  color='green', label=diff_label)
        ax2.set_xlim(min(bragg_x)-0.1, max(bragg_x)+0.1)
        ax2.legend(loc=1)
    else:
        middle_text(ax2, Bragg_text)


        
 ## ----- Partials -----
    partials_data_file = path_stem + "_PDFpartials.csv"
    ax1.set_title('Partials')
    ax1.set_xlabel(r'r($\AA$)')
    ax1.set_ylabel('g(r)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    
    if os.path.exists(partials_data_file):

        if kwargs.get("partials")=="Combine_partials":

            labels = kwargs.get("partials_labels")
            partials = kwargs.get("partials_matrix")
            
            for i in range(len(labels)-1):
                ax1.plot(partials[:,0],partials[:,i+1], label=labels[i+1])
            
            # ax1.set_xlim(right = max(partials[:,0])+0.1)
        
        else:
            with open(partials_data_file) as f:
                ncols = len(f.readline().split(','))

            partial_labels = np.genfromtxt(partials_data_file, delimiter=",", dtype='str', max_rows=1)
            partial_data = np.loadtxt(partials_data_file, delimiter=",", skiprows=1, usecols=range(0,ncols))

            x_par = partial_data[:,0]  # r(Ang)

            Va_cols = np.char.find(partial_labels, 'Va')
            non_Va_cols = np.where(Va_cols == -1)[0].tolist()

            if kwargs.get("remove_Va"):
               partial_labels = partial_labels[non_Va_cols] 
               partial_data = partial_data[:,non_Va_cols] 
            
            ncols = partial_data.shape[1]

            for n in range(1,ncols):
                PDF_par = partial_data[:,n]
                par_label = partial_labels[n]
                ax1.plot(x_par, PDF_par, label=par_label)

            # ax1.set_xlim(right = max(x_par)+0.1)
        
        ax1.set_ylim(bottom=0,)    
        ax1.legend(loc=1, ncol=3)
        


 ## ----- Chi2 -----
    # Create a list of the log files present in the directory
    l_c = 0
    log_check = stem+'-*.log'
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, log_check):
            if l_c==0: log_files = np.array([file])
            else: log_files = np.append(log_files, file)
            l_c+=1

    n_logs = len(log_files)

    log_data_File = os.path.join(path,log_files[0])

    # Set the labels
    ax3.set_xlabel('moves generated')
    ax3.set_ylabel(r'chi$^2$')
    ax3.set_title(r'moves and chi$^2$')
    # ax2_chi.set_ylabel('moves generated')
    ax2_chi.set_ylabel('moves')
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())

    if os.path.exists(log_data_File):

        # Time, moves_acc, moves_gen, fit1, fit2, fit3, etc.
        chi2_labels = np.genfromtxt(log_data_File, dtype='str', max_rows=1)
        chi2_data = np.genfromtxt(log_data_File, skip_header=2)
        # print(chi2_data.shape)

        # Set the moves and time
        if len(chi2_data.shape) == 0: pass
        elif len(chi2_data.shape) == 1:
            time_rmc = np.array([chi2_data[0]])
            move_acc = np.array([chi2_data[1]])
            move_gen = np.array([chi2_data[2]])
        else:
            time_rmc = chi2_data[:,0]
            move_acc = chi2_data[:,1]
            move_gen = chi2_data[:,2]

        moves_log = np.array([max(move_gen)])

        if n_logs > 1:

            for ln in range(1,n_logs):

                log_data_File = os.path.join(path,log_files[ln])

                # Time, moves_acc, moves_gen, fit1, fit2, fit3, etc. 
                chi2_data_log = np.genfromtxt(log_data_File, skip_header=2)

                # Set the moves and time
                if chi2_data_log.shape[0] == 0: pass
                else:
                    if len(chi2_data_log.shape) == 1:
                        # time_rmc_i = np.array([chi2_data_log[0]])
                        # move_acc_i = np.array([chi2_data_log[1]])
                        # move_gen_i = np.array([chi2_data_log[2]])

                        time_rmc_log = np.array([chi2_data_log[0]]) + max(time_rmc)
                        move_acc_log = np.array([chi2_data_log[1]]) + max(move_acc)
                        move_gen_log = np.array([chi2_data_log[2]]) + max(move_gen)
                    else:
                        # time_rmc_i = chi2_data_log[:,0]
                        # move_acc_i = chi2_data_log[:,1]
                        # move_gen_i = chi2_data_log[:,2]

                        time_rmc_log = chi2_data_log[:,0] + max(time_rmc)
                        move_acc_log = chi2_data_log[:,1] + max(move_acc)
                        move_gen_log = chi2_data_log[:,2] + max(move_gen)
                    
                    # m_t = max(time_rmc); print(m_t)
                    # m_a = max(move_acc); print(m_a)
                    # m_g = max(move_gen); print(m_g)

                    # if max(time_rmc_i) < m_t: pass
                    # elif max(time_rmc_i) > m_t: time_rmc_i -= m_t
                    
                    # if max(move_acc_i) < m_t: pass
                    # elif max(move_acc_i) > m_t: move_acc_i -= m_t
                    
                    # if max(move_gen_i) < m_t: pass
                    # elif max(move_gen_i) > m_t: move_gen_i -= m_t
                    
                    # time_rmc_log = time_rmc_i + max(time_rmc)
                    # move_acc_log = move_acc_i + max(move_acc)
                    # move_gen_log = move_gen_i + max(move_gen)

                    moves_log = np.append(moves_log,max(move_gen_log))

                    chi2_data = np.row_stack((chi2_data,chi2_data_log))

                    time_rmc = np.append(time_rmc,time_rmc_log)
                    move_acc = np.append(move_acc,move_acc_log)
                    move_gen = np.append(move_gen,move_gen_log)


        num_data = len(chi2_labels) - 3 #first 3 columns are time and moves

        # Set and plot the moves on its axis
        ax2_chi.plot(move_gen,move_acc,'r',marker='o',ms=4 , label=chi2_labels[1])
        ax2_chi.plot(move_gen,move_gen,'k',marker='o',ms=4 , label=chi2_labels[2])
        min_y = min(np.min(move_gen), np.min(move_acc)) 
        max_y = max(np.max(move_gen), np.max(move_acc))
        moves_del = 1*(max_y - min_y)
        min_y = min_y-0.001*moves_del 
        max_y = max_y+0.03*moves_del
        # ax2_chi.set_ylim(bottom=min_y,top=max_y)
        ax2_chi.set_xlim(left=min_y,right=max_y)
        ax2_chi.legend(loc=1, title='moves')
        
        # Plot each of the fit chi2 sets
        max_chi=-1000
        min_chi= 1000
        for i in range(0,num_data):
            fit_col = i + 3
            try:
                if len(chi2_data.shape) == 1:
                    chi_fit_plot = chi2_data[fit_col]
                else:
                    chi_fit_plot = chi2_data[:,fit_col]

                max_chi = max(max_chi,np.max(chi_fit_plot))
                min_chi = min(min_chi,np.min(chi_fit_plot))
            
                ax3.semilogy(move_gen,chi_fit_plot,marker='o',ms=4, label=chi2_labels[fit_col])
                # ax3.plot(move_gen,chi_fit_plot,marker='o',ms=4, label=chi2_labels[fit_col])
            except: pass
        ax3.legend(loc=2, title=r'chi$^2$')

        # Add vertical black lines corresponding to each RMC run
        moves_log_y = np.zeros(len(moves_log))
        moves_log_y[:] = 0.1
        for i in range(len(moves_log)): ax3.axvline(moves_log[i], color='k')
        # ax3.semilogy(moves_log,moves_log_y, marker='|', linewidth=0, color='k', markersize=1000)

        ax2_chi.format_coord = make_format(ax2_chi, ax3)

        moves_del_min = 0.01*(max_chi - min_chi)
        if moves_del_min<0: moves_del_min=0
        moves_del_max = 0.10*(max_chi - min_chi)
        ax3.set_ylim(bottom=min_chi-moves_del_min,top=max_chi+moves_del_max)

    from matplotlib.ticker import ScalarFormatter as SF
    f = SF(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits((0, 3))
    ax3.yaxis.set_major_formatter(f)
    ax3.xaxis.set_major_formatter(f)
    ax2_chi.yaxis.set_major_formatter(f)
    ax2_chi.xaxis.set_major_formatter(f)
            
    
 ## ----- Final touches
    cursor1 = Cursor(ax2_chi, useblit=True, color='k', linewidth=1)
    cursor2 = Cursor(ax4, useblit=True, color='k', linewidth=1)
    multi = MultiCursor(fig.canvas, (ax, ax1), useblit=True, color='k', lw=1)
        
    plt.suptitle('RMCProfile results: ' + stem, fontsize=12, fontweight='bold')
    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07, hspace=0.4)
    
    
    # return fig, cursor1, multi
    return fig, cursor1, cursor2, multi





	

