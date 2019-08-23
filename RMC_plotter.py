

# A simple and effective plotting code for assessing the fits from RMCProfile refinements


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.widgets import MultiCursor
import os

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = '11'


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def rmcplot(path, stem, pdf_type, recip_type):

    """
    Parameters
    ----
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
    """

    path_stem = path + "\\" + stem

    fig = plt.figure(figsize=cm2inch(50,25)) #width, height
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

    if os.path.exists(PDF_data_file):
        pdf_data = np.loadtxt(PDF_data_file, delimiter=",", skiprows=1) 

        pdf_r    = pdf_data[:,0]
        pdf_fit  = pdf_data[:,1]
        pdf_data = pdf_data[:,2]
        pdf_diff = pdf_data - pdf_fit
        diff_offset = abs(min(pdf_fit)) + abs(max(pdf_diff)) + 0.05
        pdf_diff = pdf_diff - diff_offset

        # Black line through difference
        ax.plot([10],[-diff_offset], '_',color='k',markersize=10000)

        ax.plot(pdf_r,pdf_data, 'bo', markerfacecolor='white', label = "PDF (Expt)")
        ax.plot(pdf_r,pdf_fit,  color='red', label='PDF (RMC)')
        ax.plot(pdf_r,pdf_diff,  color='green', label='diff')
        ax.set_xlim(min(pdf_r)-0.1, max(pdf_r)+0.1)
        ax.set_xlabel('r(A)')
        ax.set_ylabel('G(r)')
        ax.set_title('PDF fit')
        ax.legend()


    
    ## ----- Recip -----
    if recip_type=="fq":
        Recip_Data_File = path_stem + "_FQ1.csv"
    elif recip_type=="sq":
        Recip_Data_File = path_stem + "_SQ1.csv"
    
    if os.path.exists(Recip_Data_File):
        rec_data = np.loadtxt(Recip_Data_File, delimiter=",", skiprows=1) 

        rec_q    = rec_data[:,0]
        rec_fit  = rec_data[:,1]
        rec_data = rec_data[:,2]
        rec_diff = rec_data - rec_fit
        diff_offset = abs(min(rec_fit)) + abs(max(rec_diff)) + 0.05
        rec_diff = rec_diff - diff_offset

        # Black line through difference
        ax4.plot([10],[-diff_offset], '_',color='k',markersize=10000)

        ax4.plot(rec_q,rec_data, 'bo-', markerfacecolor='white', label = "Recip.(Expt)")
        ax4.plot(rec_q,rec_fit,  color='red', label='Recip.(RMC)')
        ax4.plot(rec_q,rec_diff,  color='green', label='diff')
        ax4.set_xlim(min(rec_q)-0.1, max(rec_q)+0.1)
        ax4.set_xlabel('Q(A^-1)')
        ax4.set_ylabel('F(Q),S(Q),etc.')
        ax4.set_title('Reciprocal space fit')
        ax4.legend()



    ## ----- Bragg -----
    Bragg_Data_File = path_stem + "_bragg.csv"
    if os.path.exists(Bragg_Data_File):
        bragg_data = np.loadtxt(Bragg_Data_File, delimiter=",", skiprows=1) 

        bragg_x    = bragg_data[:,0]
        bragg_exp  = bragg_data[:,1]
        bragg_fit  = bragg_data[:,2]
        bragg_diff = bragg_exp - bragg_fit
        diff_offset = abs(min(bragg_fit)) + abs(max(bragg_diff)) + 0.05
        bragg_diff = bragg_diff - diff_offset

        ax2.plot(bragg_x,bragg_exp, 'b', marker='x', markersize=4, label='Bragg (Expt)', linewidth=1.5)
        ax2.plot(bragg_x,bragg_fit,  color='red', label='Bragg (RMC)')
        ax2.plot(bragg_x,bragg_diff,  color='green', label='diff')
        ax2.set_xlim(min(bragg_x)-0.1, max(bragg_x)+0.1)
        ax2.set_xlabel('2Theta, TOF, etc.')
        ax2.set_ylabel('Intensity')
        ax2.set_title('Bragg fit')
        ax2.legend()


        
    ## ----- Partials -----
    partials_data_file = path_stem + "_PDFpartials.csv"
    if os.path.exists(partials_data_file):

        with open(partials_data_file) as f:
            ncols = len(f.readline().split(','))

        partial_labels = np.genfromtxt(partials_data_file, delimiter=",", dtype='str', max_rows=1)
        partial_data = np.loadtxt(partials_data_file, delimiter=",", skiprows=1, usecols=range(0,ncols))

        x_par = partial_data[:,0]  # r(Ang)

        for n in range(1,ncols):
            PDF_par = partial_data[:,n]
            par_label = partial_labels[n]
            ax1.plot(x_par, PDF_par, label=par_label)
        
        ax1.set_ylim(bottom=0,)
        ax1.legend(loc=0)



    ## ----- Chi2 -----
    log_num = 0
    log_num_str = str(log_num).zfill(2)
    log_data_File = path_stem + "-" + log_num_str + ".log"
        
    if os.path.exists(log_data_File):

        # Time, moves_acc, moves_gen, fit1, fit2, fit3, etc.
        chi2_labels = np.genfromtxt(log_data_File, dtype='str', max_rows=1)
        chi2_data = np.genfromtxt(log_data_File, skip_header=2)

        # Set the moves and time
        time_rmc = chi2_data[:,0]
        move_acc = chi2_data[:,1]
        move_gen = chi2_data[:,2]

        moves_log = np.array([max(move_gen)])

        log_num = 1
        while log_num < 100:
            
            log_str = str(log_num).zfill(2)

            log_data_File = path_stem + "-" + log_str + ".log"

            if os.path.exists(log_data_File):
            
            	# Time, moves_acc, moves_gen, fit1, fit2, fit3, etc. 
            	chi2_data_log = np.genfromtxt(log_data_File, skip_header=2)

            	# Set the moves and time
            	time_rmc_log = chi2_data_log[:,0] + max(time_rmc)
            	move_acc_log = chi2_data_log[:,1] + max(move_acc)
            	move_gen_log = chi2_data_log[:,2] + max(move_gen)

            	moves_log = np.append(moves_log,max(move_gen_log))

            	chi2_data = np.row_stack((chi2_data,chi2_data_log))

            	time_rmc = np.append(time_rmc,time_rmc_log)
            	move_acc = np.append(move_acc,move_acc_log)
            	move_gen = np.append(move_gen,move_gen_log)

            else:
            	break

            log_num+=1

        num_data = len(chi2_labels) - 3 #first 3 columns are time and moves

        # Set the chi2 axis 
        ax3.set_xlabel('m_generated')
        ax3.set_ylabel('chi^2')
        
        # Set and plot the moves on its axis
        ax2_chi.set_ylabel('m_accepted')
        ax2_chi.plot(move_gen,move_acc,'r' , label=chi2_labels[1])
        ax2_chi.plot(move_gen,move_gen,'k' , label=chi2_labels[2])
        ax2_chi.set_ylim(bottom=-1000,)
        ax2_chi.legend()
        
        # Plot each of the fit chi2 sets
        for i in range(0,num_data):
            fit_col = i + 3
            chi_fit_plot = chi2_data[:,fit_col]
            
            ax3.semilogy(move_gen,chi_fit_plot, label=chi2_labels[fit_col])
            ax3.legend()

        moves_log_y = np.zeros(len(moves_log))
        moves_log_y[:] = 0.1
        ax3.semilogy(moves_log,moves_log_y, marker='|', linewidth=0, color='k', markersize=1000)

        ax3.set_ylim(bottom=-1000,)		
        plt.title('chi^2')
            
    
    ## ----- Final touches
    cursor = Cursor(ax2_chi, useblit=True, color='k', linewidth=1)
    multi = MultiCursor(fig.canvas, (ax, ax1), useblit=True, color='k', lw=1)
        
    plt.suptitle('RMCProfile results: ' + stem, fontsize=12, fontweight='bold')
    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07, hspace=0.4)
    
    
    return fig, cursor, multi






	

