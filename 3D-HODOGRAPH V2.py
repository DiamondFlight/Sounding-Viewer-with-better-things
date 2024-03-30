print("loading libraries")

import os
import pip
import sharppy
import sharppy.io
from sharppy.io import *
from sharppy.io import uwyo_decoder
import numpy as np
import io
from io import StringIO
import math
import random

import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo
import matplotlib.transforms as transforms
import matplotlib
from matplotlib.axes import Axes
import matplotlib.axis as maxis
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
from matplotlib.projections import register_projection
import matplotlib.spines as mspines
from matplotlib.ticker import MultipleLocator, NullFormatter, ScalarFormatter
from matplotlib.pyplot import title
import matplotlib.patches as patches
import matplotlib.transforms as transforms
import numpy as np
import sharppy.sharptab as tab
import matplotlib as plt
import matplotlib.pyplot as plt

nameOFthing = input("Enter Location: ")

#14061619.OAX
#14061619modify.OAX
#Fiddle.txt

fName = "Fiddle.txt"

if nameOFthing == "":
    fName = "14061619.OAX"
elif nameOFthing == "mod":
    fName = "14061619modify.OAX"
elif nameOFthing == "fun":
    fName = "Fiddle.txt"
else:
    fName = "14061619.OAX"
fName = "forever sounding.txt"
    

def showSkew():
    print("load lib")
    import warnings # Silence the warnings from SHARPpy
    warnings.filterwarnings("ignore")
    import sharppy.plot.skew as skew
    from matplotlib.ticker import ScalarFormatter, MultipleLocator
    from matplotlib.collections import LineCollection
    import matplotlib.transforms as transforms
    import matplotlib.pyplot as plt
    from datetime import datetime
    import numpy as np
    from matplotlib import gridspec
    from sharppy.sharptab import winds, utils, params, thermo, interp, profile
    from sharppy.io.spc_decoder import SPCDecoder
    import matplotlib.pyplot as plt
    print("loaded lib")

    def decode(filename):

        dec = SPCDecoder(filename)

        if dec is None:
            raise IOError("Could not figure out the format of '%s'!" % filename)

        # Returns the set of profiles from the file that are from the "Profile" class.
        profs = dec.getProfiles()
        stn_id = dec.getStnId()

        for k in list(profs._profs.keys()):
            all_prof = profs._profs[k]
            dates = profs._dates
            for i in range(len(all_prof)):
                prof = all_prof[i]
                new_prof = profile.create_profile(pres=prof.pres, hght=prof.hght, tmpc=prof.tmpc, dwpc=prof.dwpc, wspd=prof.wspd, \
                                                  wdir=prof.wdir, strictQC=False, profile='convective', date=dates[i])
                return new_prof, dates[i], stn_id 
     
    FILENAME = fName

    prof, time, location = decode(FILENAME)
    sfc = prof.pres[prof.sfc]
    # Bounds of the pressure axis 
    pb_plot=1050
    pt_plot=100
    dp_plot=10
    plevs_plot = np.arange(pb_plot,pt_plot-1,-dp_plot)
    # Open up the text file with the data in columns (e.g. the sample OAX file distributed with SHARPpy)
    title = time.strftime('%Y%m%d/%H%M') + ' ' + location + '   (Observed)'
    print("located Plot")
    # Set up the figure in matplotlib.
    fig = plt.figure(figsize=(14, 10),facecolor=(0.1,0.1,0.1))
    gs = gridspec.GridSpec(4,9, width_ratios=[1,2,0.7,0.15,0.15,0.1,0.75,0.75,0.01],height_ratios=[1,1,1,1])
    ax = plt.subplot(gs[0:3, 0:3], projection='skewx', facecolor=(0,0,0))####
    ax.tick_params(labelcolor=(1,1,1,1))

    #MODIFIED
    
    skew.draw_mixing_ratio_lines(ax)
    skew.draw_dry_adiabats(ax, prof)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(prof.tmpc[np.where(~prof.tmpc.mask)[0].min()]+2, prof.pres[np.where(~prof.tmpc.mask)[0].min()] + 75,  f'{round(thermo.ctof(prof.tmpc[np.where(~prof.tmpc.mask)[0].min()]))}', verticalalignment='bottom', color=(1,0,0), fontsize=14, weight="bold")
    #ax.text(prof.tmpc[np.where(~prof.tmpc.mask)[0].min()], 0.01,  f'{round(thermo.ctof(prof.tmpc[np.where(~prof.tmpc.mask)[0].min()]))}', color=(1,0,0), transform=trans, fontsize=14, weight="bold")
    ax.text(prof.dwpc[np.where(~prof.tmpc.mask)[0].min()]+2, prof.pres[np.where(~prof.tmpc.mask)[0].min()] + 75, f'{round(thermo.ctof(prof.dwpc[np.where(~prof.tmpc.mask)[0].min()]))}', color=(0,1,0),  fontsize=14, weight="bold")
    skew.draw_title(ax, title)
    skew.draw_heights(ax, prof)
    skew.draw_effective_inflow_layer(ax,prof)
    skew.plot_sig_levels(ax, prof)

    
    ax.set_title(f'{title}',weight="bold",color=(1,1,1),loc="left")
    ax.grid(True)
    plt.grid(True)

    # Plot the background variables
    presvals = np.arange(1000, 0, -10)
    
    ax.semilogy(prof.vtmp[~prof.dwpc.mask], prof.pres[~prof.dwpc.mask], 'r--')
    ax.semilogy(prof.wetbulb[~prof.dwpc.mask], prof.pres[~prof.dwpc.mask], 'c-')
    ax.semilogy(prof.tmpc[~prof.tmpc.mask], prof.pres[~prof.tmpc.mask], 'r', lw=2)
    ax.semilogy(prof.dwpc[~prof.dwpc.mask], prof.pres[~prof.dwpc.mask],  'g', lw=2)

    
    

    # Plot the parcel trace, but this may fail.  If it does so, inform the user.
    try:
        ax.semilogy(prof.mupcl.ttrace, prof.mupcl.ptrace, 'w--')
    except:
        print("Couldn't plot parcel traces...")

    # Highlight the 0 C and -20 C isotherms.
    l = ax.axvline(0, color='b', ls='--')
    l = ax.axvline(-20, color='b', ls='--')
    #11787


    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.set_yticks(np.linspace(100,1000,10))
    ax.set_ylim(1050,100)

    # Plot the hodograph data.
    #inset_axes = skew.draw_hodo_inset(ax, prof)
    #skew.plotHodo(inset_axes, prof.hght, prof.u, prof.v, color='r')
    #inset_axes.text(srwind[0], srwind[1], 'RM', color='r', fontsize=8)
    #inset_axes.text(srwind[2], srwind[3], 'LM', color='b', fontsize=8)

    # Draw the wind barbs axis and everything that comes with it.
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlim(-50,50)
    ax2 = plt.subplot(gs[0:3,3])###
    ax3 = plt.subplot(gs[3,0:2])####
    skew.plot_wind_axes(ax2)
    skew.plot_wind_barbs(ax2, prof.pres, prof.u, prof.v)
    srwind = params.bunkers_storm_motion(prof)
    gs.update(left=0.05, bottom=0.05, top=0.95, right=1, wspace=0.025)

    # Calculate indices to be shown.  More indices can be calculated here using the tutorial and reading the params module.
    p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    
    sfc_1km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p1km)
    sfc_6km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
    srh3km = winds.helicity(prof, 0, 3000., stu = srwind[0], stv = srwind[1])
    srh1km = winds.helicity(prof, 0, 1000., stu = srwind[0], stv = srwind[1])
    scp = params.scp(prof.mupcl.bplus, prof.right_esrh[0], prof.ebwspd)
    stp_cin = params.stp_cin(prof.mlpcl.bplus, prof.right_esrh[0], prof.ebwspd, prof.mlpcl.lclhght, prof.mlpcl.bminus)
    stp_fixed = params.stp_fixed(prof.sfcpcl.bplus, prof.sfcpcl.lclhght, srh1km[0], utils.comp2vec(prof.sfc_6km_shear[0], prof.sfc_6km_shear[1])[1])
    ship = params.ship(prof)

    

    # A routine to perform the correct formatting when writing the indices out to the figure.
    def fmt(value, fmt='int'):
        if fmt == 'int':
            try:
                val = int(value)
            except:
                val = str("M")
        else:
            try:
                val = round(value,1)
            except:
                val = "M"
        return val

    

    def getHGZCape(prf):
        t,b = params.hgz(prf)
        prz = prf.pres

        #pcl = cape(prof, pres=prof.pres[i], tmpc=prof.tmpc[i], dwpc=prof.dwpc[i])
        print(t,b)

        lN = np.where((prz >= b) & (prz <= t))
        hP, hT, hD = prf.pres[lN],prf.tmpc[lN],prf.dwpc[lN]

        tCape = 0
        gTest = params.cape(prf, pbot=b, ptop=t, pres=hP, tmpc=hT, dwpc=hD)
        print(gTest,gTest.bplus)

        for pC in hP:
            tT = interp.temp(prf, pC)
            tTD = interp.dwpt(prf, pC)

            print(params.cape(prf, pres=pC, tmpc=tT, dwpc=tTD))
            tCape += params.cape(prf, pbot=t, ptop=b, pres=pC, tmpc=tT, dwpc=tTD)

        print(tCape+" CAPE")

    
    # Setting a dictionary that is a collection of all of the indices we'll be showing on the figure.
    # the dictionary includes the index name, the actual value, and the units.
    indices = {'SBCAPE': [fmt(prof.sfcpcl.bplus), 'J/kg', 1],\
               'SBCIN': [fmt(prof.sfcpcl.bminus), 'J/kg', 2],\
               'SBLCL': [fmt(prof.sfcpcl.lclhght), 'm AGL', 3],\
               'SBLFC': [fmt(prof.sfcpcl.lfchght), 'm AGL', 4],\
               'SBEL': [fmt(prof.sfcpcl.elhght), 'm AGL', 5],\
               'SBLI': [fmt(prof.sfcpcl.li5), 'C', 6],\
               'MLCAPE': [fmt(prof.mlpcl.bplus), 'J/kg', 7],\
               'MLCIN': [fmt(prof.mlpcl.bminus), 'J/kg', 8],\
               'MLLCL': [fmt(prof.mlpcl.lclhght), 'm AGL', 9],\
               'MLLFC': [fmt(prof.mlpcl.lfchght), 'm AGL', 10],\
               'MLEL': [fmt(prof.mlpcl.elhght), 'm AGL',11],\
               'MLLI': [fmt(prof.mlpcl.li5), 'C', 12],\
               'MUCAPE': [fmt(prof.mupcl.bplus), 'J/kg', 13],\
               'MUCIN': [fmt(prof.mupcl.bminus), 'J/kg', 14],\
               'MULCL': [fmt(prof.mupcl.lclhght), 'm AGL', 15],\
               'MULFC': [fmt(prof.mupcl.lfchght), 'm AGL', 16],\
               'MUEL': [fmt(prof.mupcl.elhght), 'm AGL', 17],\
               'MULI': [fmt(prof.mupcl.li5), 'C', 18],\
               '0-1 km SRH': [fmt(srh1km[0]), 'm2/s2', 19],\
               '0-1 km Shear': [fmt(utils.comp2vec(sfc_1km_shear[0], sfc_1km_shear[1])[1]), 'kts', 20],\
               '0-3 km SRH': [fmt(srh3km[0]), 'm2/s2', 21],\
               '0-6 km Shear': [fmt(utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1]), 'kts', 22],\
               'Eff. SRH': [fmt(prof.right_esrh[0]), 'm2/s2', 23],\
               'EBWD': [fmt(prof.ebwspd), 'kts', 24],\
               '3CAPE': [fmt(prof.sfcpcl.b3km), 'J/kg', 25],
               'DCAPE': [fmt(params.dcape(prof)[0]),'J/kg',26],
               '6CAPE': [fmt(prof.sfcpcl.b6km),'J/kg',27],\
               'MU MPL':[fmt(prof.mupcl.mplhght),'m AGL',28],\
               'MU FZL':[fmt(prof.mupcl.p0c),'m AGL',29],
               'PWAT': [round(prof.pwat, 2), 'inch', 30],\
               'LowRH': [round(params.mean_relh(prof, ptop=prof.mupcl.lclpres)), '%', 31],
               'MidRH': [round(params.mean_relh(prof, (prof.pres[prof.get_sfc()]-150), (prof.pres[prof.get_sfc()]-350))), '%', 32],
               'ConvT': [round(thermo.ctof(params.convective_temp(prof))), 'F', 33],\
               'MaxT': [round(thermo.ctof(params.max_temp(prof))), 'F', 34],\
               'K-index': [fmt(params.k_index(prof)), '', 35],\
               'STP(fix)': [fmt(stp_fixed, 'flt'), '', 36],\
               'SHIP': [fmt(ship, 'flt'), '', 37],\
               'SCP': [fmt(scp, 'flt'), '', 38],\
               'STP(cin)': [fmt(stp_cin, 'flt'), '', 39]}

    # List the indices within the indices dictionary on the side of the plot.
    trans = transforms.blended_transform_factory(ax.transAxes,ax.transData)
    def keyF(tick):
        for i in indices:
            if indices[i][2]==tick:
                return i
    # Write out all of the indices to the figure.
    print("##############")
    print("   INDICES    ")
    print("##############")
    string = ''
    keys = np.sort(list(indices.keys()))
    x = 0
    counter = 0
    txtColor =  (1,1,1)

    coloumns = ("CAPE","CIN","LCL","LFC","EL","LI")
    rows = ("SB","ML","MU")

    cellTxt = []
    stor = []

    for i in range(3):
        tempT = []
        tempS = []
        for m in range(6):
            key = keyF(i*6+m+1)
            tempT.append(str(indices[key][0])+' '+indices[key][1])
            tempS.append(indices[key][0])
        stor.append(tempS)
        cellTxt.append(tempT)
    
    cellTxt.reverse()
    

    tbl = plt.table(cellText=cellTxt,rowLabels=rows,cellLoc='center',colLoc='center',colLabels=coloumns,colWidths=[0.2,0.2,0.2,0.2,0.2,0.2],loc='top',edges="open",fontsize=17,bbox=[0,0.5,1,0.5])
    ax3.add_patch(patches.Rectangle(xy=(-0.3,0.5),width=1.3,height=0.5,linewidth=3,color=(0.5,0.5,0.5,0.5),fill=False))
    ax3.add_patch(patches.Rectangle(xy=(-0.1,0.89),width=1,height=0.01,linewidth=2,color=(0.5,0.5,0.5,0.5),fill=False))

    def get_color_for_Index(i,m):
        if m == 0:
            if stor[i-1][m] <= 1000:
                return (1,1,1)
            elif 1000 < stor[i-1][m] and stor[i-1][m] <= 2000:
                return (1,1,0)
            elif 2000 < stor[i-1][m] and stor[i-1][m] <= 3000:
                return (1,0.5,0)
            elif 3000 < stor[i-1][m] and stor[i-1][m] <= 4000:
                return (1,0,0)
            elif 4000 < stor[i-1][m] and stor[i-1][m] <= 5000:
                return (0.5,0,1)
            elif 5000 < stor[i-1][m] and stor[i-1][m] <= 6000:
                return (1,0,1)
            elif 6000 < stor[i-1][m] and stor[i-1][m] <= 7000:
                return (0.5,0,1)
            elif 7000 < stor[i-1][m] and stor[i-1][m] <= 8000:
                return (1,0.9,0.9)
            elif 8000 < stor[i-1][m] and stor[i-1][m] <= 9000:
                return (0.5,0.1,0.1)
            elif 9000 < stor[i-1][m] and stor[i-1][m] <= 10000:
                return (0,0.5,1)
            else:
                return (0,0.5,0)
        elif m == 1:
            if stor[i-1][m] >= -25:
                return (0,1,0)
            elif -25 > stor[i-1][m] and stor[i-1][m] >= -50:
                return (1,1,0)
            elif -50 > stor[i-1][m] and stor[i-1][m] >= -75:
                return (1,0.5,0)
            elif -75 > stor[i-1][m] and stor[i-1][m] >= -150:
                return (1,0,0)
            elif -150 > stor[i-1][m] and stor[i-1][m] >= -300:
                return (0.9,0.3,0.1)
            else:
                return (0.6,0.5,0.5)
        elif m == 4:
            if stor[i-1][m] <= 10000:
                return (1,1,1)
            elif 10000 > stor[i-1][m] and stor[i-1][m] >= 12000:
                return (0,1,0)
            elif 12000 > stor[i-1][m] and stor[i-1][m] >= 14000:
                return (0,0.5,1)
            else:
                return (0,1,1)
        else:
            return (1,1,1)

            
    
    for i in range(1,4):
        for m in range(0,6):
            colr = get_color_for_Index(i,m)
            tbl[i,m].set(facecolor=(0,0,0))
            tbl[i,m].get_text().set(x=0,y=151,color=colr,weight="bold")
        
    for i in range(0,6):
        tbl[0,i].get_text().set(color=(1,1,1),weight="bold")
    for i in range(1,4):
        tbl[i,-1].get_text().set(color=(0.9,0.8,0.7),weight="bold")
    xx=0
    for tick in range(len(indices)):
        if tick == 18:
            counter = 0
        key = keyF(tick+1)
        n = indices[key][0]
        if key == "SCP":
            if n <= 0:
                txtColor = (0,1,0)
            elif 0 < n or n <= 6:
                txtColor = (0.5,0.5,0.5)
            elif 6 < n and n <= 10:
                txtColor = (1,1,0)
            elif 10 < n and n <= 15:
                txtColor = (1,0.5,0)
            elif 15 < n and n <= 22:
                txtColor = (1,0,0)
            elif 22 < n and n <= 40:
                txtColor = (1,0,1)
            elif 40 < n and n <= 60:
                txtColor = (0.5,1,1)
            elif 60 < n and n <= 100:
                txtColor = (0,0,1)
            else:
                txtColor = (1,1,1)
        elif key == "STP(fix)" or key == "STP(cin)" or key == "SHIP":
            if n <= 0:
                txtColor = (0.5,0.8,1)
            elif 0 < n and n <= 1:
                txtColor = (0.5,0.5,0.5)
            elif 1 < n and n <= 2:
                txtColor = (1,0.5,0.5)
            elif 2 < n and n <= 3:
                txtColor = (0.8,0,0.5)
            elif 3 < n and n <= 4:
                txtColor = (1,0.5,0)
            elif 4 < n and n <= 5:
                txtColor = (1,0.5,0.5)
            elif 5 < n and n <= 6:
                txtColor = (1,0,1)
            elif 6 < n and n <= 10:
                txtColor = (0,1,0.5)
            elif 10 < n and n <= 15:
                txtColor = (1,1,0.5)
            elif 15 < n and n <= 20:
                txtColor = (0.2,0.5,1)
            elif 20 < n and n <= 40:
                txtColor = (0.6,0,0.6)
            elif 40 < n and n <= 70:
                txtColor = (1,0.2,0.2)
            elif 70 < n and n <= 100:
                txtColor = (0.75,0.25,1)
            else:
                txtColor = (1,1,1)
        else:
            txtColor = (1,1,1)
        string = string + key + ': ' + str(indices[key][0]) + ' ' + indices[key][1] + '\n'
        
        if tick <= 17:
            txtColor = (0,0,0,0)
        else:
            ax3.text(xx+0.2, 0.43-counter/7, string, verticalalignment='top', horizontalalignment='right', transform=ax3.transAxes, fontsize=8, color=txtColor, weight="bold")
        ax3.text(x, 0.5-counter/6, string, verticalalignment='top', transform=ax3.transAxes, fontsize=11, color=(0,0,0,0))
        string = ''
        txtColor = (1,1,1)
        print((key + ": " + str(indices[key][0]) + ' ' + indices[key][1]))
        if (counter < 4 and tick <= 17) or (counter < 4 and tick >= 18):
            counter += 1
            continue
        else:
            counter = 0
            #ax3.text(x, 1, string, verticalalignment='top', transform=ax3.transAxes, fontsize=11, color=txtColor)
            string = ''
            if tick > 17:
                xx += 0.225
            x += 0.3
    ax3.set_axis_off()

    # Show SARS matches (edited for Keith Sherburn)
    try:
        supercell_matches = prof.supercell_matches
        hail_matches = prof.matches 
    except:
        supercell_matches = prof.right_supercell_matches
        hail_matches = prof.right_matches

    print()
    print("#############")
    print(" SARS OUTPUT ")
    print("#############")
    for mtype, matches in zip(['Supercell', 'Hail'], [supercell_matches, hail_matches]):
        print(mtype)
        print('-----------')
        if len(matches[0]) == 0:
            print("NO QUALITY MATCHES")
        for i in range(len(matches[0])):
            try:
                print(matches[0][i] + ' ' + matches[1][i])
            except:
                print("COULDNT LOAD MATCH")
        print("Total Loose Matches:", matches[2])
        print("# of Loose Matches that met Criteria:", matches[3])
        print("SVR Probability:", matches[4])
        print() 

    # Finalize the image formatting and alignments, and save the image to the file.
    gs.tight_layout(fig)
    fn = time.strftime('%Y%m%d.%H%M') + '_' + location + '.png'
    fn = fn.replace('/', '')
    #print(("SHARPpy Skew-T image output at:", fn))
    #plt.savefig(fn, bbox_inches='tight', dpi=180)

    #@@@@@@@@@
    
    import sharppy.plot.skew as skew
    import mpl_toolkits.mplot3d.art3d as art3d

    #fig = plt.figure(figsize=(6.5875, 6.2125),facecolor=(0.2,0.2,0.2))
    #ax = fig.add_subplot(111,projection='3d',facecolor=(0.2,0.2,0.2))

    bx = plt.subplot(gs[0:2,6:9],projection='3d',facecolor=(0,0,0))
    bx.plot(4,3)
    bx.set(xlabel='U component (m/s)', ylabel='V component (m/s)', zlabel='Height (m)')
    bx.xaxis.label.set_color((1,1,1,1))
    bx.xaxis.label.set_fontsize(15)
    bx.yaxis.label.set_fontsize(15)
    bx.yaxis.label.set_color((1,1,1,1))
    bx.zaxis.label.set_fontsize(15)
    bx.zaxis.label.set_color((1,1,1,1))

    bx.tick_params(labelcolor=(1,1,1,1))

    bx.set_xlim(-60, 60)
    bx.set_ylim(-60, 60)
    bx.set_zlim(0, 12000)
    bx.set_facecolor((0,0,0))

    bx.text(0,90,0,"N",fontsize=8,color="white",horizontalalignment='center')
    bx.text(90,0,0,"E",fontsize=8,color="white",horizontalalignment='center')
    bx.text(0,-90,0,"S",fontsize=8,color="white",horizontalalignment='center')
    bx.text(-90,0,0,"W",fontsize=8,color="white",horizontalalignment='center')

    for i in range(10,90,10):
        aC = Circle((0,0),i,color=(0.8,0.8,0.8),alpha=.3,fill=False)
        bx.add_patch(aC)
        art3d.pathpatch_2d_to_3d(aC, z=0, zdir="z")
        if i%10 == 0 and i <= 50:
            bx.text(-i,2,0,str(i),color=(1,1,1,1),fontsize=8,horizontalalignment='center')
        bx.add_artist(aC)
    mot = params.bunkers_storm_motion(prof)

    aC = Circle((mot[0],mot[1]),4,color=(1,0,0),alpha=.3,fill=False)
    bx.add_patch(aC)
    art3d.pathpatch_2d_to_3d(aC, z=0, zdir="z")
    bx.text(mot[0],mot[1],5,"RM",color=(1,1,1,1),fontsize=20,horizontalalignment='center')
    bx.add_artist(aC)

    def createHodo():
        from mpl_toolkits.axes_grid.inset_locator import inset_axes
        
        imset = inset_axes(bx,width=5,height=5,loc=1)

        imset.get_xaxis().set_visible(False)
        imset.get_yaxis().set_visible(False)
        
        for i in range(10,90,10):
            am = Circle((0,0),i,color='k',alpha=.3,fill=False)
            if i%10 == 0 and i <= 50:
                imset.text(-i,2,str(i),fontsize=8,horizontalalignment='center')
            imset.add_artist(am)
            
        imset.set_xlim(-60,60)
        imset.set_ylim(-60,60)
        imset.axhline(y=0, color='k')
        imset.axvline(x=0, color='k')

        return imset


    def getAngle(point):
        return math.atan2(point[1]-mot[1],point[0]-mot[0])
    storage = []
    # Routine to plot the hodograph in segments (0-3 km,3-6 km, etc.)
    def plotHodo(axes, h, u, v, color='k'):
        for color, min_hght in zip([(1,0,0),(1,0,1),(0,1,0),(1,1,0),(0,1,1)], [[0,3000],[0,1000],[3000.1,6000],[6000.1,9000],[9000.1,12000]]):
            below_12km = np.where((h <= min_hght[1]) & (h >= min_hght[0]))[0]
            
            if len(below_12km) == 0:
                continue
            below_12km = np.append(below_12km, below_12km[-1] + 1)
            #axes.plot([mot[0]],[mot[1]],color=(1,0,1),lw=10)
            # Try except to ensure missing data doesn't cause this routine to fail.
            try:
                axes.plot(u[below_12km][~u.mask[below_12km]],v[below_12km][~v.mask[below_12km]],h[below_12km][~u.mask[below_12km]], color=color, lw=4)


                if len(storage) == 0:
                    n = np.where((h <= 12000) & (h >= 0))[0]
                    storage.append(u[n][~u.mask[n]])
                    storage.append(v[n][~v.mask[n]])
                    storage.append(h[n][~u.mask[n]])
            except:
                continue

    plotHodo(bx, prof.hght, prof.u, prof.v, color='r')

    axs = plt.subplot(gs[2,6:7],xmargin=-0.1,ymargin=-0.1,facecolor=(0.1,0.1,0.1))
    axs.plot(4,3,c=(1,1,1))

    #axs.set(xlabel='streamwise %', ylabel='height (m/s)')
    axs.set_xlim(0,1)
    axs.set_ylim(0,2.5)
    axs.tick_params(labelcolor=(1,1,1,1))

    axs.xaxis.label.set_fontsize(15)
    axs.yaxis.label.set_fontsize(15)

    #x0 = np.array([[0.8,0],[0.05,0.5],[0.02,1],[0.44,1.5],[0.75,2],[0.7,2.5]])
    #x1 = np.array([[0.2,0],[0.4,0.5],[0.8,1],[0.04,1.5],[0.05,2],[0.1,2.5]])

    #x0 = np.array([[0,0],[1,1]])
    #y0 = [0,1]

    axs.xaxis.set_tick_params(zorder=3)

    plt.rcParams["axes.axisbelow"] = True


    #axs.yaxis.label.set_color((1,1,1,1))
    plt.rcParams['axes.titley'] = 1.0    # y is in axes-relative coordinates.
    plt.rcParams['axes.titlepad'] = -14  # pad is in points...
    axs.tick_params(axis="y",direction="in", pad=-22)
    axs.tick_params(axis="x",direction="in", pad=-15)

    axs.set_title("Vorticity ($s^{-1}$)",color=(1,1,1))

    pPlot = [
        [],
        [],
        [],
        [],
        [],
        [],
        []
    ]

    vortRand = random.randint(2500,10000)/10000

    def plotVorticity(axes, h, u, v):
        hgt = h
        velo = v
        ang = u
        
        for itr in range(len(storage[1])):
            hgtMod = float(storage[2][itr])/4800
            angBn = math.degrees(getAngle([storage[0][itr],storage[1][itr]]))

            if angBn < 0:
                while angBn <0:
                    angBn+=360
            if angBn > 360:
                while angBn > 360:
                    angBn -= 360
                
            vorP = 1 - abs(135-angBn)/(225)

            #MIXING RATIO
            Mratio = interp.mixratio(prof, interp.pres(prof, storage[2][itr]))
            
            try:
                pPlot[0].append(vorP)
                pPlot[1].append(hgtMod)
                pPlot[2].append(vorP-random.randint(0,100)/705)
                pPlot[3].append(hgtMod)
                pPlot[4].append(Mratio)
            except:
                continue



    plotVorticity(axs,prof.hght,prof.u,prof.v)
    axs.plot(pPlot[2], pPlot[3], label='ω',zorder=2,color=(0.5,0.5,0.5))
    axs.plot(pPlot[0], pPlot[1], label='ωₛ',color=(1,0,0),zorder=1)
    locs, labels = plt.xticks()  # Get the current locations and labels.
    plt.xticks(np.arange(0, 1, step=0.2),[0,0.2,0.4,0.6,0.8],zorder=3)  # Set label locations.
    #plt.xticks(np.arange(3), ['Tom', 'Dick', 'Sue'])  # Set text labels.
    #plt.xticks([0, 1, 2], ['January', 'February', 'March'], rotation=20)  # Set text labels and properties.



    axs.legend(loc="lower right", fontsize="7")


    axe = plt.subplot(gs[2,7:8],xmargin=-0.1,ymargin=-0.1,facecolor=(0.1,0.1,0.1))
    axe.plot(4,3,c=(1,1,1))

    #axs.set(xlabel='streamwise %', ylabel='height (m/s)')
    axe.set_xlim(0,20)
    axe.set_ylim(0,2.5)
    axe.tick_params(labelcolor=(1,1,1,1))

    axe.xaxis.label.set_fontsize(15)
    axe.yaxis.label.set_fontsize(15)

    #x0 = np.array([[0.8,0],[0.05,0.5],[0.02,1],[0.44,1.5],[0.75,2],[0.7,2.5]])
    #x1 = np.array([[0.2,0],[0.4,0.5],[0.8,1],[0.04,1.5],[0.05,2],[0.1,2.5]])

    #x0 = np.array([[0,0],[1,1]])
    #y0 = [0,1]

    axe.xaxis.set_tick_params(zorder=3)

    plt.rcParams["axes.axisbelow"] = True


    #axs.yaxis.label.set_color((1,1,1,1))
    plt.rcParams['axes.titley'] = 1.0    # y is in axes-relative coordinates.
    plt.rcParams['axes.titlepad'] = -14  # pad is in points...
    axe.tick_params(axis="y",direction="in", pad=-22)
    axe.tick_params(axis="x",direction="in", pad=-15)
    plt.xticks(np.arange(0, 20, step=5),[5,10,15,20],zorder=3)

    axe.set_title("Mix Ratio (g/kg)",color=(1,1,1))
    axe.plot(pPlot[4], pPlot[3], label='w',color=(0,1,0))
    axe.text(5, 1.8, f'Eff RH: {round(params.mean_relh(prof,prof.ebottom,prof.etop))}%', color=(1,1,1), fontsize=10)
    axe.legend(loc="lower right", fontsize="7")
    #Theta E
    thetaE = prof.get_thetae_profile()
    idx = np.where( prof.pres > 400 )[0]
    
    tEmin = prof.thetae[idx].min()-10
    tEmax = prof.thetae[idx].max()+10
    print(tEmin,tEmax)
    ace = plt.subplot(gs[3,6:7],xmargin=-0.1,ymargin=-0.1,facecolor=(0.1,0.1,0.1))
    ace.plot(1,1,c=(1,1,1))
    ace.set_xlim(tEmin,tEmax)
    ace.set_ylim(0,600)
    mk1 = prof.thetae.mask
    mk2 = prof.pres.mask
    mk = np.maximum(mk1,mk2)
    pTE = prof.pres[~mk]
    tTE = prof.thetae[~mk]

    for i in range(pTE.shape[0] -1):
        if pTE[i] > 400:
            pPlot[5].append(tTE[i])
            pPlot[6].append(1000 - pTE[i])

    #axs.set(xlabel='streamwise %', ylabel='height (m/s)')
    ace.tick_params(labelcolor=(1,1,1,1))

    ace.xaxis.label.set_fontsize(15)
    ace.yaxis.label.set_fontsize(15)
    ace.plot(pPlot[5],pPlot[6], color=(1,0,0))

    #x0 = np.array([[0.8,0],[0.05,0.5],[0.02,1],[0.44,1.5],[0.75,2],[0.7,2.5]])
    #x1 = np.array([[0.2,0],[0.4,0.5],[0.8,1],[0.04,1.5],[0.05,2],[0.1,2.5]])

    #x0 = np.array([[0,0],[1,1]])
    #y0 = [0,1]

    ace.xaxis.set_tick_params(zorder=3)

    plt.rcParams["axes.axisbelow"] = True


    #axs.yaxis.label.set_color((1,1,1,1))
    plt.rcParams['axes.titley'] = 1.0    # y is in axes-relative coordinates.
    plt.rcParams['axes.titlepad'] = -14  # pad is in points...
    ace.tick_params(axis="y",direction="in", pad=-22)
    ace.tick_params(axis="x",direction="in", pad=-15)
    plt.xticks(np.arange(round(tEmin/10)*10 + 20, round(tEmax/10)*10, step=10),zorder=3)
    plt.yticks(np.arange(100, 500, step=100),[900,800,700,600])

    ace.set_title("Theta-E vs Pressure",color=(1,1,1))

    ax4 = plt.subplot(gs[3,2:4])####
    ax4.set_axis_off()
    ax4.add_patch(patches.Rectangle(xy=(0,0.5),width=1,height=0.5,linewidth=2,color=(0.5,0.5,0.5,0.5),fill=False))
    
    lRates = [
        ['SFC-3km AGL',round(10*params.lapse_rate(prof, 0, 3000, False))/10],
        ['3-6km AGL',round(10*params.lapse_rate(prof, 3000, 6000, False))/10],
        ['850-500mb',round(10*params.lapse_rate(prof, 850, 500, True))/10],
        ['700-500mb',round(10*params.lapse_rate(prof, 700, 500, True))/10],
    ]

    for i in range(len(lRates)):
        txtColor = (1,1,1)

        if lRates[i][1] < 5:
            txtColor = (0.8,0.8,0.8)
        elif lRates[i][1] > 5 and lRates[i][1] <= 6.8:
            txtColor = (0.7,1,0.7)
        elif lRates[i][1] > 6.8 and lRates[i][1] <= 8:
            txtColor = (0.9,0.9,0)
        elif lRates[i][1] > 8 and lRates[i][1] <= 9:
            txtColor = (1,0,0)
        elif lRates[i][1] > 9 and lRates[i][1] <= 9.8:
            txtColor = (0.9,0.4,0.7)
        elif lRates[i][1] > 9.8:
            txtColor = (1,1,1)
        ax4.text(0.05, 0.95-i*0.125, f'{lRates[i][0]} LR: {lRates[i][1]} C/km', verticalalignment='top', transform=ax4.transAxes, fontsize=8, color=txtColor, weight="bold")
    

spc_file = open(fName, 'r').read()


Wyf = open('input.txt','r').read()



print("loaded")

def parseSPC(spc_file):
    """
        This function will read a SPC-style formatted observed sounding file,
        similar to that of the 14061619.OAX file included in the SHARPpy distribution.

        It will return the pressure, height, temperature, dewpoint, wind direction and wind speed data
        from that file.
    """
    ## read in the file
    newDat = spc_file.split('\n')

    def findType(fnd):
        xi = -1
        for x in newDat:
            xi+=1
            if fnd in x:
                return xi
    titleDex = findType("%TITLE%")
    startDex = findType("%RAW%")+1
    endDex = findType("%END%")

    pltTitle = newDat[titleDex+1].split()
    ##

    fullDat = '\n'.join(newDat[startDex : endDex][:])

    soundDat = StringIO(fullDat)

    p, h, T, Td, wdir, wspd = np.genfromtxt( soundDat, delimiter=',', comments='%', unpack=True)


    Date = pltTitle
    
    return p,h,T,Td,wdir,wspd,Date


pres, hght, tmpc, dwpc, wdir, wspd, Date = parseSPC(spc_file)

def formatDate(date):
    fStr = ""
    fStr = fStr + f"{date[0]}    "
    fStr = fStr + f"{date[1][2:4]}"
    fStr = fStr + f"/{date[1][4:6]}"
    fStr = fStr + f"/20{date[1][0:2]}"
    fStr = fStr + f"     {date[1][7:9]}:{date[1][9:]}z      (OBSERVED)"
    return fStr

graphDate = formatDate(Date)


prof = profile.create_profile(profile='default', pres=pres, hght=hght, tmpc=tmpc, \
                                    dwpc=dwpc, wspd=wspd, wdir=wdir, missing=-9999, strictQC=False)



showSkew()





plt.show()
