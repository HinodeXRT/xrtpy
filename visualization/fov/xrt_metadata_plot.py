
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import matplotlib.colors as mcolors


def get_pcol(ndata):
    cmap = plt.get_cmap('rainbow')
    ndim = ndata#+1
    col_arr = np.linspace(0.0,1.0,ndim)
    col_lis = cmap(col_arr)#[:-1])
    return col_lis

def lab_locs0(xc, yc):
    xlocs = np.asarray([np.min(xc.value), np.max(xc.value), np.max(xc.value), np.min(xc.value)])
    ylocs = np.asarray([np.min(yc.value), np.min(yc.value), np.max(yc.value), np.max(yc.value)])
    hal   = ['left'    , 'right'   , 'right'   , 'left'   ]
    val   = ['bottom'  ,'bottom'   , 'top'     , 'top'     ]
    return xlocs, ylocs, hal, val

def mk_label_locs(xlis0, ylis0, dycheck = 90.0):
    xlis, ylis, ylis0 = np.asarray(xlis0), np.asarray(ylis0), np.asarray(ylis0)
    xpos = np.min(xlis) + 300.0
    dist_bool = True
    ymaxi = np.argmax(ylis)
    for j in range(100):#while dist_bool:
        good_bool = True
        if False:#good_bool:
            print()
            print(ylis)
        for i in range(len(ylis)):
            if (i != ymaxi):
                dyy = np.abs(ylis - ylis[i])
                dyy[i] = 10.0**3.0
                if (np.min(dyy) < dycheck):
                    #print(dyy)
                    ylis[i] = ylis[i] - dycheck
                    #ylisb[i] = ylisb[i] - dycheck*2.0
                    good_bool = False
                    break
            
        if good_bool:
            #print('yo!')
            dist_bool = False
        if not dist_bool:
            break
    dyl = ylis - ylis0
    ylisb = ylis0 + 0.5*dyl
    return np.ones(len(xlis))*xpos, ylis, ylisb

def fix_col(col_in, dcol_v):
    colv = mcolors.to_rgb(col_in)
    colv2 = np.asarray([colv[0], colv[1],colv[2]])
    sum_v = dcol_v - np.sum(colv2)
    dcol = 1.0-colv2
    #plt.fill_between([i,i+1.0],[0.0]*2,[1.0]*2,color=colv2)
    colv3 = colv2 + dcol*(sum_v)/np.sum(dcol)
    return colv3        




def interactive_plot(filter_selection, 
                     zoom_v, 
                     tzoom, 
                     time_ind, 
                     mini_tline,
                     inputs):
    

    # interactive plot function

    xca = inputs['xca']
    yca     = inputs['yca'] 
    rsun_p  = inputs['rsun_p'] 
    zla     = inputs['zla'] 
    lon_arr = inputs['lon_arr'] 
    yla     = inputs['yla'] 
    xmeta   = inputs['xmeta'] 
    flen    = inputs['flen'] 
    t0d     = inputs['t0d'] 
    trp     = inputs['trp'] 
    dt_bool = inputs['dt_bool'] 
    fov_lis = inputs['fov_lis']
    fill_bool  = inputs['fill_bool'] 
    filter_lis = inputs['filter_lis'] 
    time_lis0  = inputs['time_lis0']
    col_vals   = inputs['col_vals']
    time_lis_abs = inputs['time_lis_abs']

    plt.close('all')

    #Set up Figure

    if mini_tline:
        fig, axs2 = plt.subplots(3,1,figsize=(8,11), gridspec_kw={'height_ratios': [16, 1, 4]})
        axs = [axs2[0],axs2[2],axs2[1]]
        jlen = 2
    else:
        fig, axs = plt.subplots(2,1,figsize=(8,10), gridspec_kw={'height_ratios': [16, 4]})
        jlen = 1
        
    ##plot the solar limb and grid lines
    axs[0].plot(xca, yca,'w-')
    axs[0].plot(np.asarray([-1.0,1.0])*rsun_p,np.zeros(2),'w-',lw=0.5)
    for j in range(len(lon_arr)):
        axs[0].plot(xca*lon_arr[j],yca,'w-',lw=0.5)
        axs[0].plot([-yla[j],yla[j]],[zla[j]]*2,'w-',lw=0.5)
        axs[0].plot([-yla[j],yla[j]],[-zla[j]]*2,'w-',lw=0.5)
        

    # Zoom Button Affects this
    time_arr = np.linspace(xmeta.metadata['TIME_RANGE']*tzoom[0],xmeta.metadata['TIME_RANGE']*tzoom[1],100)

    # In case we want to change the center point
    #x0r, y0r = x0*rsun_p, y0*rsun_p

    
    alpha_arr = np.ones(flen)*0.1
    alpha_arr2 = np.zeros(flen)
    
    if filter_selection == 'All':
        alpha_arr = np.ones(flen)
        fbool = True
    else:
        fbool = False
        fi = filter_lis.index(filter_selection)
        alpha_arr[fi] = 1.0
        alpha_arr2[fi] = 0.3

    # Current time for timeline
    time_ii = (np.timedelta64(int(np.round(time_arr[time_ind])),'s')+t0d)

    
    for j in range(jlen):
        axs[1+j].plot([time_ii]*2, [-0.5, flen-0.5], 'w-',lw=4,zorder=1000)
        axs[1+j].plot([time_ii]*2, [-0.5, flen-0.5], 'b-',zorder=1001)  
        
    xll, yll = [], []
    dtl = []
    fplis = []
    alpha_plis = []
    col_pvals = []
    # Range for the zoomed in timeline
    trp2 = [np.timedelta64(int(np.round(time_arr[0])),'s')+t0d,(np.timedelta64(int(np.round(time_arr[-1])),'s')+t0d)]

    if mini_tline:
        axs[2].plot([trp2[0]]*2, [-0.5, flen-0.5], 'w-',lw=2,zorder=200)  
        axs[2].plot([trp2[1]]*2, [-0.5, flen-0.5], 'w-',lw=2,zorder=200)  
    
        if fill_bool:
            axs[2].fill_between([trp[0],trp2[0]],[-0.5]*2,[flen-0.5]*2,color='k',alpha = 0.5,zorder=100)
            axs[2].fill_between([trp[1],trp2[1]],[-0.5]*2,[flen-0.5]*2,color='k',alpha = 0.5,zorder=100)
    
        lout = axs[2].plot([trp[0],trp2[0]],[-8.0,-0.5],'w-',lw=0.5)
        lout2 = axs[2].plot([trp[1],trp2[1]],[-8.0,-0.5],'w-',lw=0.5)

    
        lout[0].set_clip_on(False)
        lout2[0].set_clip_on(False)

    #print(fov_lis[0][0])
    
    for i in range(flen):
        for j in range(jlen):
            axs[1+j].plot(trp,i*np.ones(2),'w-',alpha = 1.0,lw=0.5,zorder=1)
            
        #Find nearest obs
        ti = np.argmin(np.abs(np.asarray(time_lis0[i])-time_arr[time_ind]))

        

        #Getting the locs the frame
        xc, yc = fov_lis[i][ti][0][0],fov_lis[i][ti][0][1] #fov_lis[filter_i][time_i][0][x,y]'
        dtv = time_lis0[i][ti]-time_arr[time_ind]
        # Only if the frame is visible
        if (np.abs(dtv) < 60.0):
            label_bool = False
            if fbool:#(i == fi):
                label_bool = True
            elif (i == fi):
                label_bool = True
            else:
                label_bool = False
                
            if label_bool:
                #this silly function just gives you the four corners of frame, 
                #plus the horizontal and vertical alignment, if you wanted to put the labels inside the frame
                xlocs, ylocs, hal, val = lab_locs0(xc, yc)# this function gives
        
                #get distance between image and current time
                dtl.append(dtv)
                # corners [BL, BR, TR, TL]
                xll.append(xlocs[-2]) #Top-right corner
                yll.append(ylocs[-2])
                fplis.append(filter_lis[i])
                alpha_plis.append(alpha_arr[i])
                col_pvals.append(col_vals[i])
        #time list for the image points for each filter
        tpl = np.array(time_lis_abs[i],dtype='datetime64[ms]')#[]


        #plot the image times on the time line(s)
        
        for j in range(jlen):
            axs[1+j].plot(tpl,np.ones(len(time_lis0[i]))*i, 'h',alpha = alpha_arr[i],color=col_vals[i])    

        # Only plot the frame if the observation is within 60s of the current time
        if (np.abs(dtv) < 60.0):
            axs[0].plot(xc, yc,'--',alpha = alpha_arr[i],color=col_vals[i],lw=1)
            if fill_bool:
                axs[0].fill_between([xc[0].value,xc[1].value],[yc[0].value]*2,[yc[2].value]*2,color=col_vals[i],alpha = alpha_arr2[i])
        
            #highlight the current frame in time line(only if it is plotted)
            axs[1].plot(tpl[ti], i, 'h',alpha = alpha_arr[i],color='white',markersize=15)
            axs[1].plot(tpl[ti], i, 'h',alpha = alpha_arr[i],color=col_vals[i],markersize=10)

        #Add label to the time line for each filter
        axs[1].text(trp2[-1]+np.timedelta64(10,'s'), i, filter_lis[i],alpha = alpha_arr[i],
                fontsize=12,color=col_vals[i],horizontalalignment='left',verticalalignment='center',fontweight='bold')
   
    #This function makes sure the filter labels on the plot don't overlap and are legible  

    if (len(dtl) > 0):
        xlis, ylis, ylisb = mk_label_locs(xll, yll)
    
    for i in range(len(dtl)):
        
        if (np.abs(dtl[i]) < 60.0):
             #Pre-amble handles the delta T for each frame
            dyp = 50.0
            col_rgb = col_pvals[i]#fix_col(col_vals[i], 2.0)
            if (dtl[i] >= 0.0):
                sstr = '+'
            else:
                sstr = '-'
            if dt_bool:
                col_tim = [np.clip(np.abs(dtl[i])/20.0,0.0,1.0),0.0,0.0]
                col_tim = [1.0,1.0-np.clip(np.abs(dtl[i])/20.0,0.0,1.0),1.0-np.clip(np.abs(dtl[i])/20.0,0.0,1.0)]
                dtstr = sstr + str(np.round(10.0*np.abs(dtl[i]))/10.0)
    
                
       
            axs[0].plot([xlis[i], xll[i]], [ylis[i]+dyp,ylisb[i]],'-',alpha = alpha_plis[i],color=col_rgb,lw=1)
            if dt_bool:
                axs[0].text( xlis[i], ylis[i]+dyp, '('+dtstr+'s)',alpha = alpha_plis[i],
                        fontsize=10,color=col_tim,horizontalalignment='right',verticalalignment='bottom',fontweight='bold')
        
        
        
            axs[0].text( xlis[i], ylis[i]+dyp, fplis[i],alpha = alpha_plis[i],#*0.6,
                        fontsize=15,color=col_rgb,horizontalalignment='left',fontweight='bold',verticalalignment='bottom')
    for j in range(jlen):        
        axs[1+j].set_yticks([])

    plot_col = 'dimgrey'
    for j in range(jlen+1):   
        axs[j].set_facecolor(plot_col)
    
    fig.patch.set_facecolor(plot_col)
    axs[0].set_title(str(time_ii),fontsize=20,color='white')
    axs[1].set_xlim(trp2)
    if (jlen > 1):
        axs[2].set_xlim(trp)
    
    locator = mdates.AutoDateLocator(minticks=4, maxticks=9)
    formatter = mdates.ConciseDateFormatter(locator)
    axs[1].xaxis.set_major_locator(locator)
    axs[1].xaxis.set_major_formatter(formatter)
    if (jlen > 1):
        locator = mdates.AutoDateLocator(minticks=4, maxticks=9)
        formatter = mdates.ConciseDateFormatter(locator)
        axs[2].xaxis.set_major_locator(locator)
        axs[2].xaxis.set_major_formatter(formatter)
        axs[2].set_ylim([-0.75, flen-0.25])
    axs[0].set_xlim(np.asarray([-1300.0,1300.0])/zoom_v)#+x0r)
    axs[0].set_ylim(np.asarray([-1300.0,1300.0])/zoom_v)#+y0r)
    
    #axs[2].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    for ax in axs:
        ax.tick_params(color='white', labelcolor='white')
        for spine in ax.spines.values():
            spine.set_edgecolor('white')

    plt.show()#zebra
    #plt.close('all')

