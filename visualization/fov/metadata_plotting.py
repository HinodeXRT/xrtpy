
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
cmap_b = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#01EA79',
                                                                  '#01EAFF',
                                                                  '#952EFC',
                                                                  '#F51A8C',
                                                                  '#FF9138'])

cmap_w = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#019488',
                                                                  '#003066',
                                                                  '#51389B',#'#3D318B',
                                                                  '#B42D7E',
                                                                  '#FF504D',
                                                                  '#E59500'])

def get_pcol(ndata, n_mode = True):
    cmap = cmap_b#plt.get_cmap('rainbow')
    cmap2 = cmap_w
    ndim = ndata#+1
    col_arr = np.linspace(0.0,1.0,ndim)
    col_lis, col_lis2 = cmap(col_arr), cmap2(col_arr)#[:-1])
    return col_lis, col_lis2

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
                     night_mode,
                     inputs,
                     nfilt,
                     alt_bool):
    #print('new interact2')
    plt.close('all')
    #alt_bool = False
    if mini_tline:
        if alt_bool:
            fig, axs2 = plt.subplots(3,1,figsize=(9,8+nfilt*0.4), gridspec_kw={'height_ratios': [16, 1, nfilt]})
            axs = [axs2[0],axs2[2],axs2[1]]
        else:
            fig, axs2 = plt.subplot_mosaic([['t1','t1','im'],['t1','t1','im'],['t2','t2','im']], figsize=(18,6))#,['im','t0','t0'],['im','t3','t3']]
            #axs2['t0'].axis('off')
            #axs2['t3'].axis('off')
            axs = [axs2['im'],axs2['t1'],axs2['t2']]
        jlen = 2
        
    else:
        if alt_bool:
            fig, axs = plt.subplots(2,1,figsize=(9,7+nfilt*0.4), gridspec_kw={'height_ratios': [16, nfilt]})

        else:
            fig, axs2 = plt.subplot_mosaic([['t1','t1','im'],['t1','t1','im']], figsize=(18,6))#,['im','t0','t0'],['im','t3','t3']]
            #axs2['t0'].axis('off')
            #axs2['t3'].axis('off')
            #axs2['t2'].axis('off')
            axs = [axs2['im'],axs2['t1']]
        jlen = 1

    axs, fig = fov_plotter(axs, fig,
                filter_selection, 
                zoom_v, 
                tzoom, 
                time_ind, 
                mini_tline,
                night_mode,
                jlen,
                inputs,
                alt_bool)
    plt.show()


def get_frames(filter_selection, 
                     zoom_v, 
                     tzoom, 
                     time_ind, 
                     mini_tline,
                     night_mode,
                     inputs,
                     nfilt,
                     alt_bool):
    #print('new interact2')
    plt.close('all')
    if mini_tline:
        if alt_bool:
            fig, axs2 = plt.subplots(3,1,figsize=(9,8+nfilt*0.4), gridspec_kw={'height_ratios': [16, 1, nfilt]})
            axs = [axs2[0],axs2[2],axs2[1]]
        else:
            fig, axs2 = plt.subplot_mosaic([['t1','t1','im'],['t1','t1','im'],['t2','t2','im']], figsize=(18,6))#,['im','t0','t0'],['im','t3','t3']]
            #axs2['t0'].axis('off')
            #axs2['t3'].axis('off')
            axs = [axs2['im'],axs2['t1'],axs2['t2']]
        jlen = 2
    else:
        fig, axs = plt.subplots(2,1,figsize=(9,10), gridspec_kw={'height_ratios': [16, 4]})
        jlen = 1

    axs, fig = fov_plotter(axs, fig,
                filter_selection, 
                zoom_v, 
                tzoom, 
                time_ind, 
                mini_tline,
                night_mode,
                jlen,
                inputs,
                alt_bool)
    return axs, fig

def make_animation2(filter_lis_b,inputs, dmode, alt_bool, nfilt, demo, tlen = 20):
    
    axs_lis = []
    
    ti = np.linspace(0,99,tlen,dtype=int)

    for i in range(tlen):
        #print(ti[i])
        axs, fig = get_frames(filter_lis_b[0],
                        1.0,
                        [0.0,1.0],
                        ti[i],
                        demo,
                        dmode,
                        inputs,
                        nfilt,
                        alt_bool
                        )
        axs_lis.append(fig)
        plt.close()
    Xlis = []
    for i in range(tlen):
        fig = axs_lis[i]
        canvas = FigureCanvasAgg(fig)
    
    
        # Retrieve a view on the renderer buffer
        canvas.draw()
        buf = canvas.buffer_rgba()
            # convert to a NumPy array
        X = np.asarray(buf)
        Xlis.append(X)
        plt.close('all')
    if alt_bool:
        fig, ax = plt.subplots(figsize=(8,11))
    else:
        fig, ax = plt.subplots(figsize=(18,6))
    plt.axis('off')
    #plt.axis("tight")  # gets rid of white border
    #plt.axis("image")
    ims = []
    for i in range(tlen):
        # make a Figure and attach it to a canvas.
        X = Xlis[i]
        #print(i)
        im = ax.imshow(X, animated=True)
        #plt.axis("tight")
        if i == 0:
            ax.imshow(X)  # show an initial one first
                #plt.close()
        ims.append([im])
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)#, wspace=None, hspace=None)
    ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True,
                                repeat_delay=1000)
    return ani #plt.show()

def make_animation(filter_lis_b,inputs, dmode, alt_bool,nfilt, demo, tlen = 20):
    
    axs_lis = []
    
    ti = np.linspace(0,99,tlen,dtype=int)
    j_zoom = 1
    if demo:
        j_zoom = 4
    jlis = [[0.0,1.0],
            [0.2,0.5],
            [0.4,0.7],
            [0.1,0.2]]

        
    for j in range(j_zoom):
        for i in range(tlen):
            #print(ti[i])
            axs, fig = get_frames(filter_lis_b[0],
                            1.0,
                            jlis[j],
                            ti[i],
                            True,
                            dmode,
                            inputs,
                            nfilt,
                            alt_bool
                            )
            axs_lis.append(fig)
            plt.close()
    Xlis = []
    for i in range(len(axs_lis)):
        fig = axs_lis[i]
        canvas = FigureCanvasAgg(fig)
    
    
        # Retrieve a view on the renderer buffer
        canvas.draw()
        buf = canvas.buffer_rgba()
            # convert to a NumPy array
        X = np.asarray(buf)
        Xlis.append(X)
        plt.close('all')
    if alt_bool:
        fig, ax = plt.subplots(figsize=(8,11))
    else:
        fig, ax = plt.subplots(figsize=(18,6))
    plt.axis('off')
    #plt.axis("tight")  # gets rid of white border
    #plt.axis("image")
    ims = []
    for i in range(len(axs_lis)):
        # make a Figure and attach it to a canvas.
        X = Xlis[i]
        #print(i)
        im = ax.imshow(X, animated=True)
        #plt.axis("tight")
        if i == 0:
            ax.imshow(X)  # show an initial one first
                #plt.close()
        ims.append([im])
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)#, wspace=None, hspace=None)
    ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True,
                                repeat_delay=1000)
    return ani #plt.show()

def fov_plotter(axs, fig,
                filter_selection, 
                     zoom_v, 
                     tzoom, 
                     time_ind, 
                     mini_tline,
                     night_mode,
                     jlen,
                     inputs,
                     alt_bool):#,
                     #focus_bool):
 
    focus_bool = False#True
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
    col_vals1   = inputs['col_vals']
    col_vals2   = inputs['col_vals2']
    time_lis_abs = inputs['time_lis_abs']


    
    if night_mode:
        col_vals = col_vals1.copy()

        plot_col2 = '#0F121F'#'#272A35'#'#282C3B'#'dimgrey'
        plot_col1 = '#272A35'#'#0F121F'#'#161820'
        line_col = '#03FFFF'
        gcol = 'w'
        lab_col = 'white'
    else:
        col_vals = col_vals2.copy()
        plot_col2 = '#FFFFFF'#'#272A35'#'#282C3B'#'dimgrey'
        plot_col1 = '#E9E9E9'#'#0F121F'#'#161820'
        line_col = '#272A35'
        gcol = 'k'
        lab_col = 'black'
    
    
    
        
    ##plot the solar limb and grid lines
    axs[0].plot(xca, yca,gcol+'-')
    axs[0].plot(np.asarray([-1.0,1.0])*rsun_p,np.zeros(2),gcol+'-',lw=0.5)
    for j in range(len(lon_arr)):
        axs[0].plot(xca*lon_arr[j],yca,gcol+'-',lw=0.5)
        axs[0].plot([-yla[j],yla[j]],[zla[j]]*2,gcol+'-',lw=0.5)
        axs[0].plot([-yla[j],yla[j]],[-zla[j]]*2,gcol+'-',lw=0.5)
        

    # Zoom Button Affects this
    time_arr = np.linspace(xmeta.metadata['TIME_RANGE']*tzoom[0],xmeta.metadata['TIME_RANGE']*tzoom[1],100)

    # In case we want to change the center point
    #x0r, y0r = x0*rsun_p, y0*rsun_p

    
    alpha_arr = np.ones(flen)*0.1
    alpha_arr2 = np.zeros(flen)
    
    if filter_selection == 'All':
        fig_title = 'All Filters'
        mini_title = ''
        mini_col = lab_col
        alpha_arr = np.ones(flen)
        fbool = True
    else:
        fig_title = 'Filter :'#filter_selection
        mini_title = filter_selection
        
        fbool = False
        fi = filter_lis.index(filter_selection)
        mini_col = col_vals[fi]
        alpha_arr[fi] = 1.0
        alpha_arr2[fi] = 0.3

    # Current time for timeline
    time_ii = (np.timedelta64(int(np.round(time_arr[time_ind])),'s')+t0d)

    dtmax = 60.0

    dt_step = time_arr[1]-time_arr[0]
    if (dt_step > dtmax):
        dtmax = dt_step
    if (dt_step > 999.0):
        dtmax = 999.0
    
    for j in range(jlen):
        axs[1+j].plot([time_ii]*2, [-0.5, flen-0.5], color=line_col,lw=4,zorder=1000)
        axs[1+j].plot([time_ii]*2, [-0.5, flen-0.5], gcol+'-',lw=1,zorder=1001)
        
        #axs[1+j].plot([time_ii]*2, [-0.5, flen-0.5], 'b-',zorder=1001)  
        
    xll, yll, all, all2 = [], [], [], []
    xmm, ymm = [], []
    xll2 = []
    dtl = []
    fplis = []
    alpha_plis = []
    col_pvals = []
    # Range for the zoomed in timeline
    trp2 = [np.timedelta64(int(np.round(time_arr[0])),'s')+t0d,(np.timedelta64(int(np.round(time_arr[-1])),'s')+t0d)]
    trp3 = [np.timedelta64(int(np.round(time_arr[1])),'s')+t0d,(np.timedelta64(int(np.round(time_arr[-2])),'s')+t0d)]
    if mini_tline:
        x_vals = [trp3[0],trp2[0],trp2[0],trp3[0]]
        y_vals = [-0.6,-0.4, flen-0.6, flen-0.4]
        axs[2].plot(x_vals, y_vals,color=line_col,lw=2,zorder=200)  

        x_vals = [trp3[1],trp2[1],trp2[1],trp3[1]]
        y_vals = [-0.5,-0.5, flen-0.5, flen-0.5]
        axs[2].plot(x_vals, y_vals,color=line_col,lw=2,zorder=200)  
    
        if fill_bool:
            axs[2].fill_between([trp[0],trp2[0]],[-0.5]*2,[flen-0.5]*2,color=plot_col2,alpha = 0.5,zorder=100)
            axs[2].fill_between([trp[1],trp2[1]],[-0.5]*2,[flen-0.5]*2,color=plot_col2,alpha = 0.5,zorder=100)
            axs[2].fill_between([trp2[0],trp2[1]],[-0.5]*2,[flen-0.5]*2,color=plot_col1,alpha = 1.0,zorder=0)
        if False:
            lout = axs[2].plot([trp[0],trp2[0]],[-8.0,-0.5],'w-',lw=0.5)
            lout2 = axs[2].plot([trp[1],trp2[1]],[-8.0,-0.5],'w-',lw=0.5)

    
            lout[0].set_clip_on(False)
            lout2[0].set_clip_on(False)

    #print(fov_lis[0][0])
    minx_lis = []
    miny_lis = []
    maxx_lis = []
    maxy_lis = []
    #tobs_lis = []
    #xc_lis, yc_lis = [], []
    for i in range(flen):
        for j in range(jlen):
            axs[1+j].plot(trp,i*np.ones(2),gcol+'-',alpha = 1.0,lw=0.5,zorder=1)
            
        #Find nearest obs
        ti = np.argmin(np.abs(np.asarray(time_lis0[i])-time_arr[time_ind]))

        
        #min_x = 10.0**6.0
        #max_x = -10.0**6.0

        #min_y = 10.0**6.0
        #max_y = -10.0**6.0

        #Getting the locs the frame
        xc, yc = fov_lis[i][ti][0][0],fov_lis[i][ti][0][1] #fov_lis[filter_i][time_i][0][x,y]'
        dtv = time_lis0[i][ti]-time_arr[time_ind]
        # Only if the frame is visible
        if (np.abs(dtv) < dtmax):
            #tobs = (np.timedelta64(int(np.round(time_arr[0])),'s')+t0d)
            
            #tobs_lis.append(tobs)
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
                xmm.append(np.mean(xlocs))
                ymm.append(np.mean(ylocs))
                #get distance between image and current time
                dtl.append(dtv)
                # corners [BL, BR, TR, TL]
                if (np.mean(xlocs) > 0.0):
                    xll.append(xlocs[3]) #Top-right corner
                    yll.append(ylocs[3])
                    all.append('right')
                    all2.append('left')
                    xll2.append(xlocs[0]-300.0)
                else:
                    xll.append(xlocs[2]) 
                    yll.append(ylocs[2])
                    all.append('left')
                    all2.append('right')
                    xll2.append(xlocs[1]+300.0)
                fplis.append(filter_lis[i])
                alpha_plis.append(alpha_arr[i])
                col_pvals.append(col_vals[i])
        #time list for the image points for each filter
        tpl = np.array(time_lis_abs[i],dtype='datetime64[ms]')#[]


        #plot the image times on the time line(s)
        dttt = trp2[-1] -  (np.timedelta64(int(np.round(time_arr[-2])),'s')+t0d)
        for j in range(jlen):
            axs[1+j].plot(tpl,np.ones(len(time_lis0[i]))*i, 'h',alpha = alpha_arr[i],color=col_vals[i])    

        # Only plot the frame if the observation is within 60s of the current time
        if (np.abs(dtv) < dtmax):
            minx_lis.append(np.min(xc.value))
            miny_lis.append(np.min(yc.value))
            maxx_lis.append(np.max(xc.value))
            maxy_lis.append(np.max(yc.value))
            #xc_lis.append(np.argmin(np.abs(xc)))
            #yc_lis.append(yc)
            axs[0].plot(xc, yc,'--',alpha = alpha_arr[i],color=col_vals[i],lw=1)
            if fill_bool:
                axs[0].fill_between([xc[0].value,xc[1].value],[yc[0].value]*2,[yc[2].value]*2,color=col_vals[i],alpha = alpha_arr2[i])
        
            #highlight the current frame in time line(only if it is plotted)

            axs[1].plot(tpl[ti], i, 'h',alpha = alpha_arr[i],color='white',markersize=15)
            axs[1].plot(tpl[ti], i, 'h',alpha = alpha_arr[i],color=col_vals[i],markersize=10)

        #Add label to the time line for each filter
        #+np.timedelta64(10,'s')
        #time_ii = (np.timedelta64(int(np.round(time_arr[time_ind])),'s')+t0d)
        
        if alt_bool:
            ggg = axs[1].text(trp2[-1]+dttt, i, filter_lis[i],alpha = alpha_arr[i],
                fontsize=15,color=col_vals[i],horizontalalignment='left',verticalalignment='center',fontweight='bold')
        else:
            ggg = axs[1].text(trp2[0]-dttt, i, filter_lis[i],alpha = alpha_arr[i],
                fontsize=15,color=col_vals[i],horizontalalignment='right',verticalalignment='center',fontweight='bold')
        ggg.set_clip_on(False)
   
    #This function makes sure the filter labels on the plot don't overlap and are legible  

    if (len(dtl) > 0):
        #xlis, ylis, ylisb = [], [], []
        #for h in range(len(dtl)):
        #    xcoor = np.asarray([])
        xlis, ylis, ylisb = mk_label_locs(xll, yll)

    xlim_out = np.asarray([-1300.0,1300.0])/zoom_v
    ylim_out = np.asarray([-1300.0,1300.0])/zoom_v

    
    if (len(dtl) > 0):
        if True:#(zoom_v > 1.0):
            zoom_val = 10.0**np.abs(np.log10(zoom_v))
            xmean, ymean = np.mean(xmm), np.mean(ymm)
            intx = np.linspace(1.0,2.0,10,endpoint=True)
            int_func = np.linspace(0.0,1.0,10,endpoint=True)**0.5
            #if (zoom_v <= 1.0):
            #    int_func = np.linspace(0.0,1.0,10,endpoint=True)**0.25

            xcen = np.interp(zoom_val,intx, int_func)*xmean
            ycen = np.interp(zoom_val,intx, int_func)*ymean
            xlim_out = xlim_out + xcen 
            ylim_out = ylim_out + ycen

    if False:#focus_bool:
        if (len(dtl) > 0):
            over_p = 1.1
            xmin, xmax, ymin, ymax = np.min(np.asarray(minx_lis)),np.max(np.asarray(maxx_lis)),np.min(np.asarray(miny_lis)),np.max(np.asarray(maxy_lis))

            mxp, dxp = (xmax+xmin)/2.0, (xmax-xmin)/2.0
            myp, dyp = (ymax+ymin)/2.0, (ymax-ymin)/2.0
            if (dxp > dyp):
                dpp = dxp*over_p
            else:
                dpp = dyp*over_p
            xlim_out = mxp + (np.asarray([-dpp, dpp])/zoom_v)
            ylim_out = myp + (np.asarray([-dpp, dpp])/zoom_v)

    dyl = (ylim_out[1]-ylim_out[0])*0.33
    yh_lis = np.linspace(ylim_out[0]+dyl,ylim_out[1],flen+2)
    #print(yh_lis)
    for i in range(len(dtl)):
        
        if (np.abs(dtl[i]) < dtmax):
            tobs = (np.timedelta64(int(np.round(dtl[i])),'s')+time_ii)
             #Pre-amble handles the delta T for each frame
            dyp = 50.0
            col_rgb = col_pvals[i]#fix_col(col_vals[i], 2.0)
            if (dtl[i] >= 0.0):
                sstr = '+'
            else:
                sstr = '-'
            if dt_bool:
                col_f = np.clip(np.abs(dtl[i])*3.0/dtmax,0.0,1.0)
                col_tim = [col_f,0.0,0.0]
                if night_mode:
                    col_tim = [1.0,1.0-col_f,1.0-col_f]
                dtstr2 = str(tobs)
                dtstr2 = dtstr2[11:-4]
                dtstr = sstr + str(np.round(10.0*np.abs(dtl[i]))/10.0)
    
                
            #dxt = -2000.0
            
            if alt_bool:
                dxxx = 0.0#-2000.0
            else:
                dxxx= 0.0#-2000.0
            if True:
                axs[0].plot([xll2[i], xll[i]], [ylis[i]+dyp,ylisb[i]],'-',alpha = alpha_plis[i],color=col_rgb,lw=1)
                if dt_bool:
                    g1 = axs[0].text( xll2[i], ylis[i]+dyp, '('+dtstr+'s)',alpha = alpha_plis[i],
                            fontsize=10,color=col_tim,horizontalalignment=all2[i],verticalalignment='bottom',fontweight='bold')
                    g1.set_clip_on(True)
        
        
        
                g2 = axs[0].text( xll2[i], ylis[i]+dyp, fplis[i],alpha = alpha_plis[i],#*0.6,
                        fontsize=15,color=col_rgb, horizontalalignment=all[i],fontweight='bold',verticalalignment='bottom')
                g2.set_clip_on(True)
            else:
                axs[0].plot([xlim_out[-1], xll[i]], [yh_lis[i+1],ylisb[i]],'-',alpha = alpha_plis[i],color=col_rgb,lw=2)
    for j in range(jlen):        
        axs[1+j].set_yticks([])

    
    for j in range(jlen+1):   
        axs[j].set_facecolor(plot_col1)
    
    fig.patch.set_facecolor(plot_col2)
    #axs[0].set_title('UTC '+str(time_ii),fontsize=20,color=mini_col)
    axs[0].set_title(fig_title+mini_title , fontsize=35,color=mini_col,weight='bold')
    axs[1].set_title(''+str(time_ii),fontsize=20,color=line_col)#,loc =  'left')
    
    #axs[2].set_title('Mini Timeline',fontsize=10,color='white',loc =  'left')
    

    #fig.suptitle(fig_title+mini_title , fontsize=35,color=mini_col, x = 0.4,y = 0.98,weight='bold')#, ha = 'left'x = 0.15,
    axs[1].set_xlim(trp2)
    if (jlen > 1):
        axs[2].set_xlim(trp)
        axs[2].set_title('Timeline Range Selection',fontsize=12,color=lab_col,zorder=2000)
    
    locator = mdates.AutoDateLocator(minticks=4, maxticks=9)
    formatter = mdates.ConciseDateFormatter(locator)
    axs[1].xaxis.set_major_locator(locator)
    axs[1].xaxis.set_major_formatter(formatter)
    if (jlen > 1):
        #locator = mdates.AutoDateLocator(minticks=4, maxticks=9)
        #formatter = mdates.ConciseDateFormatter(locator)
        #axs[2].xaxis.set_major_locator(locator)
        #axs[2].xaxis.set_major_formatter(formatter)
        #axs[2].set_xticks([])
        axs[2].set_ylim([-0.75, flen-0.25])
        axs[2].axis('off')

    
    if True:
        axs[0].text( xlim_out[0]*0.95,ylim_out[-1]*0.98,'UTC '+str(time_ii),alpha = 0.7,#*0.6,
                        fontsize=15,color=mini_col, horizontalalignment='left',fontweight='bold',verticalalignment='top')
    #if xlim_out[-1] < 1500:
    #    xlim_out[-1] = 1500

    axs[0].set_xlim(xlim_out)#+x0r)
    axs[0].set_ylim(ylim_out)#+y0r)
    
    #axs[2].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    for ax in axs:
        ax.tick_params(color=lab_col, labelcolor=lab_col)#,direction='out',length=5)
        #ax.tick_params(color=lab_col, labelcolor=lab_col,direction='in',which='minor')
        for spine in ax.spines.values():
            spine.set_edgecolor(line_col)#'white')
    axs[0].xaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0].grid(alpha=0.5,ls=':')
    plt.tight_layout()
    return axs, fig

